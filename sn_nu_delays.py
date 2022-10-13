
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy.coordinates import SkyCoord
from datetime import date, datetime, timedelta

import numpy as np
import pandas as pd
import plotly.graph_objs as go
from format_json import get_json_data

from dash import Dash, dcc, html, Input, Output

EarthRadius = 6.3781e6 # in units m
Earthradkpc = EarthRadius*u.m.to(u.kpc)
isoformat = "%Y-%m-%dT%H:%M:%S.%f"

# radius = EarthRadius # in units m
uu, vv = np.mgrid[0:2*np.pi:200j, 0:np.pi:100j]
xE = EarthRadius * np.cos(uu)*np.sin(vv)
yE = EarthRadius * np.sin(uu)*np.sin(vv)
zE = EarthRadius * np.cos(vv)

# class from SNEWPDAG
class Detector:
    def __init__(self, name, lon, lat, height, sigma, bias):
        self.name = name
        self.lon = lon  # degrees
        self.lat = lat  # degrees
        self.height = height  # [m]
        self.sigma = sigma  # time resolution [s]
        self.bias = bias  # time bias [s], observed - true
        delta = np.radians(90.0 - lat)
        alpha = np.radians(lon)
        self.r = np.array([np.sin(delta) * np.cos(alpha),
                           np.sin(delta) * np.sin(alpha),
                           np.cos(delta)])
        self.loc = EarthLocation(lon=lon, lat=lat)

    def get_gcrs(self, obstime):
        # k = EarthLocation.of_site('keck')
        t = Time(obstime, format='isot')  # make sure it's in astropy Time form
        g = self.loc.get_gcrs(obstime=t)  # g.ra and g.dec
        return g

    def get_xyz(self, obstime, scale_rad=1):
        g = self.get_gcrs(obstime)
        radius = np.sqrt(self.loc.x ** 2 + self.loc.y ** 2 + self.loc.z ** 2)
        radius *= scale_rad
        codelta = np.radians(g.dec)
        alpha = np.radians(g.ra)
        sphi = np.sin(alpha)
        cphi = np.cos(alpha)
        ctheta = np.sin(codelta)
        stheta = np.cos(codelta)
        return radius * np.array([stheta * cphi, stheta * sphi, ctheta])

def get_xyz_from_gcrs(gcrs_coor, radius=10, scale_rad=1):
    """ radius in kpc
        returns x,y,z in GCRS coor, in meters
    """
    radius = (radius * u.kpc).to(u.m).value * scale_rad
    delta = gcrs_coor.dec.rad
    phi = gcrs_coor.ra.rad
    # SNx = -np.cos(phi) * np.cos(delta)
    # SNy = -np.sin(phi) * np.cos(delta)
    # SNz = -np.sin(delta)
    SNx = np.cos(phi) * np.cos(delta)
    SNy = np.sin(phi) * np.cos(delta)
    SNz = np.sin(delta)
    return radius * np.array([SNx, SNy, SNz])

#------------------------------------------------ load assets ----------------------------------------------------------
## Detectors data
detectors_df = pd.read_csv('./assets/detector_locations.csv', index_col='Unnamed: 0')

## Stars data
stars_radec_list = np.loadtxt('./assets/star_candidates.txt', usecols=(1, 2), delimiter=',')
star_random_radii = np.random.uniform(low=1.5, high=12, size=len(stars_radec_list))
candid_names = np.loadtxt('./assets/star_candidates.txt', usecols=0, delimiter=',', dtype=str)
star_names_list = []
for _name in candid_names:
    stripped_name = _name.strip()
    if stripped_name not in star_names_list:
        star_names_list.append(stripped_name)
    else:
        star_names_list.append(stripped_name + "-2")

_stars = np.column_stack((star_names_list, stars_radec_list, star_random_radii*Earthradkpc))
stars_df = pd.DataFrame(_stars, columns=('Star Name', 'RA', 'DEC', 'Rand. Dist'))
stars_df = stars_df.astype({"RA": float, "DEC": float, 'Rand. Dist':float})
stars_df.set_index('Star Name', inplace=True)

# to be displayed on hover, fixed
star_text = [name +
             f"<br>RA:{ra:.2f}" +
             f"<br>DEC:{dec:.2f}" for name, ra, dec in zip(stars_df.index, stars_df['RA'], stars_df['DEC'])]
#-----------------------------------------------------------------------------------------------------------------------

def get_detectors_at_time(obstime_str="2023-06-14T12:00:00.000000", detectors=detectors_df):
    obstime = Time(obstime_str, format='isot', scale='utc')
    df_new = detectors.copy()
    df_new['ObsTime'] = np.repeat(obstime, len(detectors))
    df_new['x (m)'], df_new['y (m)'], df_new['z (m)'] = 0, 0, 0
    # xyz_detectors = np.zeros((len(selected_detectors),3))
    for i, name in enumerate(detectors.index):
        det = detectors.loc[name]
        detobj = Detector(name=name, lon=det['Longitude (deg)'], lat=det['Latitude (deg)'],
                          height=det['Height (m)'], sigma=det['Time Uncertainty (s)'],
                          bias=det['Bias (s)'])
        df_new.loc[name, 'x (m)'], df_new.loc[name, 'y (m)'], df_new.loc[name, 'z (m)'] = detobj.get_xyz(obstime).value
        df_new.loc[name, 'RA (deg)'] = detobj.get_gcrs(obstime).ra.value
        df_new.loc[name, 'DEC (deg)'] = detobj.get_gcrs(obstime).dec.value
        df_new.loc[name, 'r (m)'] = detobj.get_gcrs(obstime).distance.value
    return df_new

def get_stars_at_time(obstime_str="2023-06-14T12:00:00.000000"):
    obstime = Time(obstime_str, format='isot', scale='utc')
    skycoors = SkyCoord(stars_df['RA'], stars_df['DEC'], frame="fk5", obstime=obstime, unit='deg')
    radii_stars = stars_df['Rand. Dist'].values

    star_xyz = get_xyz_from_gcrs(skycoors.gcrs, radius=radii_stars)
    star_dict = {'Star Name': stars_df.index,
                 'RA (deg)': stars_df['RA'],
                 'DEC (deg)': stars_df['DEC'],
                 'x (m)': star_xyz[0],
                 'y (m)': star_xyz[1],
                 'z (m)': star_xyz[2]}

    stars_xyz = pd.DataFrame.from_dict(star_dict).set_index('Star Name')
    return stars_xyz

def star_loc_at_time(star_name, obstime_str="2023-06-14T12:00:00.000000",):
    stars_xyz = get_stars_at_time(obstime_str)
    star = stars_xyz[stars_xyz.index==star_name]
    return star

####### Make dash app

### The wireframe for Earth
lines = []
line_marker = dict(color='cyan', width=3) # '#0066FF'
for i, j, k in zip(xE, yE, zE):
    lines.append(go.Scatter3d(x=i, y=j, z=k, mode='lines', line=line_marker, hoverinfo='none'))

layout = go.Layout(
    title=dict(text ='SN Detectors & Candidate Stars',
               font =dict(family='Sherif',
               size=14,
               color = 'cyan')),
    width=800,
    height=800,
    paper_bgcolor='rgba(0,0,0,1)',
    plot_bgcolor='rgba(0,0,0,1)',
    showlegend = False,
    scene=dict(
        xaxis=dict(
            gridcolor='rgb(0,0,0,0)',
            zerolinecolor='rgb(0,0,0,0)',
            showbackground=False,
            backgroundcolor='rgb(0,0,0,0)',
            visible=False,
            linecolor='rgba(0,0,0,0)'
        ),
        yaxis=dict(
            gridcolor='rgb(0,0,0,0)',
            zerolinecolor='rgb(0,0,0,0)',
            showbackground=False,
            backgroundcolor='rgb(0,0,0,0)',
            visible=False,
            linecolor='rgba(0,0,0,0)'
        ),
        zaxis=dict(
            gridcolor='rgb(0,0,0,0)',
            zerolinecolor='rgb(0,0,0,0)',
            showbackground=False,
            backgroundcolor='rgb(0,0,0,0)',
            visible=False
        ),
        # camera=dict(eye=dict(x=1.15, y=1.15, z=1.15))
    ),
    uirevision=None,
)


def scatter_selected_star(detectors, stars, exploding):
    candid = stars.loc[exploding]
    radec = f"RA:{candid['RA (deg)']}, DEC:{candid['DEC (deg)']}"
    modified_part =[ go.Scatter3d(x=[candid['x (m)']], y=[candid['y (m)']], z=[candid['z (m)']],
                                  mode='markers', hoverinfo='text',
                                  marker=dict(size=35, color='yellow', opacity=0.4),
                                  text=f"{exploding} <br>{radec}<br> EXPLODED!")]
    lines_to_detectors = []
    x_star = candid['x (m)']
    y_star = candid['y (m)']
    z_star = candid['z (m)']
    for detector in detectors.index:
        det = detectors.loc[detector]
        xs = np.array([x_star, det['x (m)']])
        ys = np.array([y_star, det['y (m)']])
        zs = np.array([z_star, det['z (m)']])
        lines_to_detectors.append(go.Scatter3d(x=xs, y=ys,z=zs, mode='lines', line = dict(color="cyan")))

    return lines_to_detectors + modified_part

def get_scatter_data(detectors, stars):
    detector_text = [name +
                      f"<br>lon: {detectors.loc[name, 'Longitude (deg)']:.2f}deg<br>" \
                      f"lat: {detectors.loc[name, 'Latitude (deg)']:.2f}deg<br>" \
                      f"Height{detectors.loc[name, 'Height (m)']:.2f}m<br>" \
                      f"RA: {detectors.loc[name, 'RA (deg)']:.2f}deg<br>" \
                      f"DEC: {detectors.loc[name, 'DEC (deg)']:.2f}deg<br>" \
                      f"r: {detectors.loc[name, 'r (m)']:.3f}m"
                      for name in detectors.index]

    det_scat =[go.Scatter3d(x=detectors['x (m)'], y=detectors['y (m)'], z=detectors['z (m)'],
                            marker=dict(size=7, symbol='x'), text=detector_text, hoverinfo='text', mode ='markers')]
    star_scat = [go.Scatter3d(x=stars['x (m)'], y=stars['y (m)'], z=stars['z (m)'], mode ='markers',
                              hoverinfo='text', text=star_text, marker=dict(size=5, color='white'))]
    return star_scat + det_scat

def plot_spheres(detectors, stars, exploding):
    globe = lines
    stars_and_detectors = get_scatter_data(detectors, stars)
    exploding_star = scatter_selected_star(detectors, stars, exploding)
    data = globe + stars_and_detectors + exploding_star
    fig = go.Figure(data=data, layout=layout)
    return fig


def delays_for_time_star(star_xyz, detectors_df, obstime_str="2023-06-14T12:00:00.000000"):
    """ For a given time, and given Star return detector delays
    """
    delayslist = []
    for detector in detectors_df.index:
        detector_xyz = detectors_df.loc[detector, ['x (m)', 'y (m)', 'z (m)']]
        distance = np.linalg.norm(star_xyz - detector_xyz)  # meters
        delayslist.append(distance / 2.9979e8)

    minval = min(np.array(delayslist))
    delays = np.array(delayslist) - minval

    timeslist = []
    for detector, delays_sec in zip(detectors_df.index, delays):
        date_obj = datetime.strptime(obstime_str, isoformat) + timedelta(seconds=delays_sec)
        timeslist.append(datetime.strftime(date_obj, isoformat))
    timeslist = np.array(timeslist)

    sortedargs = np.argsort(delays)
    delays_arr = delays[sortedargs]
    times_arr = timeslist[sortedargs]
    detectors_arr = detectors_df.index[sortedargs]
    delays_df = pd.DataFrame({'Detectors':detectors_arr, 'delays':np.round(delays_arr,6), 'times':times_arr}, index=detectors_arr,)
    delays_df = delays_df[delays_df.columns.tolist()[::-1]]
    return delays_df


def generate_table(dataframe, max_rows=20):
    ## convert df to iso time formatted
    return html.Table([
        html.Thead(
            html.Tr([html.Th(col) for col in dataframe.columns])
        ),
        html.Tbody([
            html.Tr([
                html.Td(dataframe.iloc[i][col]) for col in dataframe.columns
            ]) for i in range(min(len(dataframe), max_rows))
        ])
    ])


star_names_list = sorted(star_names_list)
kv_pairs = [{"label":name, "value":val} for val, name in enumerate(star_names_list)]
vk_pairs = {val:name for val, name in enumerate(star_names_list)}

### Make a dashboard
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = Dash(__name__, external_stylesheets=external_stylesheets)
server = app.server

explanation_text =  f"\n\nThe app calculates the relative time delays between selected_detectors \n" \
                    f"Therefore, the radial distance of a candidate star is randomly distributed \n" \
                    f"The plane wave from a star is brought to earth surface, and all delays are referenced\n" \
                    f"from the first detector interaction.\n" \
                    "see the code at [GitHub](https://github.com/KaraMelih/supernova-triangulation)"

app.layout = html.Div([
    html.Div(
        children=[
            html.H1("Supernova Arrival Time Delays", style={'text-align': 'center', 'color':'white'}),
            dcc.Markdown(children=explanation_text, style={'color':'white', 'font_size': '6px'}),
            dcc.Dropdown(id="candid_selected", options=kv_pairs, multi=False, value=123,
                         placeholder="Select a Star to Explode", style={'width': "90%"}),
            dcc.Checklist(id="Detectors", options=detectors_df.index, value=['SK', 'IC', 'KL', 'XENONnT'],
                          style={'color':'white', 'font_size': '10px'}, inline=True),
            dcc.DatePickerSingle(id='my-date-picker-single', min_date_allowed=date(2022, 1, 1),
                                 max_date_allowed=date(2030, 12, 12), initial_visible_month=date(2023, 6, 14),
                                 date=date(2023, 6, 14)),
            html.Div(id='output-container-date-picker-single'),
            # download button
            html.A('Download Data',
                    id='download-link',
                    download="rawdata.json",
                    href="",
                    target="_blank", style={'font_size': '15px'}),
            ],
        style={'display': 'inline-block', 'width': "30%", 'vertical-align': 'top', 'margin-left': '3vw',
               'margin-top': '3vw'}
    ),

    html.Div(
        children=[
            html.Div([html.H4(children='Delays in sec', style={'color':'white'}),
                      html.Table(id='delay_df', style={"font-weight": "bold", "color": "white"})])
        ],
        style={'display': 'inline-block', 'vertical-align': 'top', 'margin-left': '3vw', 'margin-top': '3vw'}
    ),
    html.Br(),
    html.Div(children=[dcc.Graph(id='my_3dscat', figure={}),],
            style={'display': 'inline-block', 'vertical-align': 'top', 'margin-left': '3vw', 'margin-top': '3vw'}),],
    style={'display': 'flex', 'backgroundColor': 'black'})


# ------------------------------------------------------------------------------
# Connect the Plotly graphs with Dash Components

@app.callback(
    [Output(component_id='my_3dscat', component_property='figure'),
     Output(component_id='delay_df', component_property='children'),
     Output('download-link', 'href'),],
    [Input(component_id='my-date-picker-single', component_property='date'),
     Input(component_id='Detectors', component_property='value'),
     Input(component_id='candid_selected', component_property='value'),])
def update_positions(obstime, selected_detectors, candid_selected):
    candid_name = vk_pairs[candid_selected]
    if obstime is not None:
        date_object = date.fromisoformat(obstime)
        date_string = date_object.strftime('%Y-%m-%d')+" 12:00:00.000000"
        date_string = datetime.strptime(date_string, "%Y-%m-%d %H:%M:%S.%f").strftime(isoformat)
    else:
        date_string = "2023-06-14T12:00:00.000000"

    detector_pos = get_detectors_at_time(date_string)
    detectors_selected = detector_pos[detector_pos.index.isin(selected_detectors)]
    star_pos = get_stars_at_time(date_string)
    fig = plot_spheres(detectors_selected, star_pos, candid_name)
    # delay df
    star_xyz = star_pos.loc[candid_name, ['x (m)', 'y (m)', 'z (m)']]
    delay_df = delays_for_time_star(star_xyz, detectors_selected, obstime_str=date_string)

    href = get_json_data(date_string, detectors_selected, delay_df)
    return fig, generate_table(delay_df), href


# ------------------------------------------------------------------------------
if __name__ == '__main__':
    app.run_server(debug=True)

