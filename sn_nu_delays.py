
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy.coordinates import SkyCoord
from datetime import date

import numpy as np
import pandas as pd
import plotly.graph_objs as go

from dash import Dash, dcc, html, Input, Output

EarthRadius = 6.3781e6 # in units m
Earthradkpc = EarthRadius*u.m.to(u.kpc)

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
        t = Time(obstime)  # make sure it's in astropy Time form
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


def append_at_obs(timestring, dataframe):
    """ append x,y,z values at a given obs time
    """
    dataframe['ObsTime'] = timestring
    ObsTime = Time(timestring, format='iso', scale='utc')
    dataframe['x (m)'] = -999
    dataframe['y (m)'] = -999
    dataframe['z (m)'] = -999
    for _name in detectors_df.index:
        detector = detectors_df.loc[_name]
        detector_obj = Detector(name=_name,
                                lon=detector['Longitude (deg)'],
                                lat=detector['Latitude (deg)'],
                                height=detector['Height (m)'],
                                sigma=detector['Time Uncertainty (s)'],
                                bias=detector['Bias (s)'])
        x, y, z = detector_obj.get_xyz(ObsTime)
        detectors_df.loc[_name, 'x (m)'] = x.value
        detectors_df.loc[_name, 'y (m)'] = y.value
        detectors_df.loc[_name, 'z (m)'] = z.value


def get_xyz_from_gcrs(gcrs_coor, radius=10, scale_rad=1):
    """ radius in kpc
        returns x,y,z in GCRS coor, in meters
    """
    radius = (radius * u.kpc).to(u.m).value * scale_rad
    delta = gcrs_coor.dec.rad
    phi = gcrs_coor.ra.rad
    SNx = -np.cos(phi) * np.cos(delta)
    SNy = -np.sin(phi) * np.cos(delta)
    SNz = -np.sin(delta)
    return radius * np.array([SNx, SNy, SNz])

#------------------------------------------------ load assets ----------------------------------------------------------
stars_radec_list = np.loadtxt('./assets/star_candidates.txt', usecols=(1, 2), delimiter=',')
candid_names = np.loadtxt('./assets/star_candidates.txt', usecols=0, delimiter=',', dtype=str)
star_names_list = []
for _name in candid_names:
    stripped_name = _name.strip()
    if stripped_name not in star_names_list:
        star_names_list.append(stripped_name)
    else:
        star_names_list.append(stripped_name + "-2")

detectors_df = pd.read_csv('./assets/detector_locations.csv', index_col='Unnamed: 0')

def get_detectors_at_time(obstime_str="2022-06-14 20:00:00.100", detectors=detectors_df):
    obstime = Time(obstime_str, format='iso', scale='utc')
    df_new = detectors.copy()
    df_new['ObsTime'] = np.repeat(obstime, len(detectors))
    df_new['x (m)'], df_new['y (m)'], df_new['z (m)'] = 0, 0, 0
    xyz_detectors = np.zeros((len(detectors),3))
    for i, name in enumerate(detectors.index):
        det = detectors.loc[name]
        detobj = Detector(name=name, lon=det['Longitude (deg)'], lat=det['Latitude (deg)'],
                          height=det['Height (m)'], sigma=det['Time Uncertainty (s)'],
                          bias=det['Bias (s)'])
        df_new.loc[name, 'x (m)'], df_new.loc[name, 'y (m)'], df_new.loc[name, 'z (m)'] = detobj.get_xyz(obstime).value
    return df_new

def get_stars_at_time(obstime_str="2022-06-14 20:00:00.100", stars_radec=stars_radec_list,
                      radii='random', star_names=star_names_list):
    assert len(star_names_list) == len(stars_radec), "length of Names and Coordinates don't match!"
    random_radii = np.random.uniform(low=1.5, high=12, size=len(star_names))
    ra = stars_radec[:, 0]
    dec = stars_radec[:, 1]
    obstime = Time(obstime_str, format='iso', scale='utc')
    skycoors = SkyCoord(ra, dec, frame="fk5", obstime=obstime, unit='deg')
    if type(radii)==str and radii=='random':
        radii_stars = Earthradkpc * random_radii
    else:
        radii_stars = np.repeat(radii, len(ra))
    star_xyz = get_xyz_from_gcrs(skycoors.gcrs, radius=radii_stars)
    star_dict = {'name':star_names, 'RA (deg)':ra, 'DEC (deg)':dec,
                 'x (m)':star_xyz[0], 'y (m)':star_xyz[1], 'z (m)':star_xyz[2]}
    star_df = pd.DataFrame.from_dict(star_dict).set_index('name')
    return star_df


def get_delays_for_time(obstime_str="2022-06-14 20:00:00.100",
                        detectors=detectors_df,
                        stars_radec=stars_radec_list):
    """ For a given time, return detector delays
        to all stars, as a dictionary contains; my_dict[star_name] => delays
    """

    stars_df = get_stars_at_time(obstime_str=obstime_str, stars_radec=stars_radec, star_names=star_names_list)
    detectors_df = get_detectors_at_time(obstime_str="2022-06-14 20:00:00.100", detectors=detectors)
    #     sn_candidate = stars_df.loc[star_name]

    title = f'{obstime_str} SN neutrino arrival delay [s]'

    delays_dict = {}
    for star in stars_df.index:
        delays_dict[star] = {}
        star_xyz = stars_df.loc[star, ['x (m)', 'y (m)', 'z (m)']]
        for detector in detectors_df.index:
            detector_xyz = detectors_df.loc[detector, ['x (m)', 'y (m)', 'z (m)']]
            distance = np.linalg.norm(star_xyz - detector_xyz)  # meters
            delays_dict[star][detector] = distance / 2.9979e8  # seconds

        delays_dict[star] = pd.DataFrame.from_dict(delays_dict[star], orient='index', columns=[title]).sort_values(
            by=title)
    return delays_dict


def get_delays_for_star(star_name,
                        delays_dict=None,
                        obstime_str="2022-06-14 20:00:00.100",
                        detectors=detectors_df,
                        stars_radec=stars_radec_list):
    if delays_dict is None:
        delays_dict = get_delays_for_time(obstime_str=obstime_str,
                                          detectors=detectors,
                                          stars_radec=stars_radec)
        star_df = delays_dict[star_name]
    else:
        star_df = delays_dict[star_name]
    return star_df


####### Make dash app

lines = []
line_marker = dict(color='cyan', width=3) # '#0066FF'
for i, j, k in zip(xE, yE, zE):
    lines.append(go.Scatter3d(x=i, y=j, z=k, mode='lines', line=line_marker, hoverinfo='none'))

layout = go.Layout(
    title='SN Detectors & Candidate Stars',
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
        camera=dict(eye=dict(x=1.15, y=1.15, z=1.15))
    ),

)

detector_text =[name + f"<br>lon: {detectors_df.loc[name, 'Longitude (deg)']:.2f}deg<br>lat: " \
                       f"{detectors_df.loc[name, 'Latitude (deg)']:.2f}deg<br>" \
                       f"Height{detectors_df.loc[name, 'Height (m)']:.2f}m" for name in detectors_df.index]
star_text = [name + f"<br>RA:{ra:.2f}" + f"<br>DEC:{dec:.2f}" for name, ra, dec in zip(star_names_list, stars_radec_list[:, 0], stars_radec_list[:, 1])]

def get_scatter_data(obstime):
    """ for a given time, get all the scatter data
    """
    detectors = get_detectors_at_time(obstime_str=obstime)
    det_scat =[go.Scatter3d(x=detectors['x (m)'], y=detectors['y (m)'], z=detectors['z (m)'],
                            marker=dict(size=7, symbol='x'), text=detector_text, hoverinfo='text', mode ='markers')]
    stars = get_stars_at_time(obstime_str=obstime)
    star_scat = [go.Scatter3d(x=stars['x (m)'], y=stars['y (m)'], z=stars['z (m)'], mode ='markers',
                              hoverinfo='text', text=star_text, marker=dict(size=5, color='white'))]
#     delays_dict = get_delays_for_time(obstime)
    return star_scat + det_scat


def plot_spheres(obstime, candid_name=None):
    data = get_scatter_data(obstime)
    if candid_name is not None:
        stars = get_stars_at_time(obstime_str=obstime)
        candid = stars.loc[candid_name]
        modified_part =[go.Scatter3d(x=[candid['x (m)']], y=[candid['y (m)']], z=[candid['z (m)']],
                                     mode='markers', hoverinfo='text',
                                     marker=dict(size=35, color='yellow', opacity=0.4),
                                     text=f"{candid_name} <br> EXPLODED!")]
        data = lines + data + modified_part
    fig = go.Figure(data=data, layout=layout)
    return fig

def generate_table(dataframe, max_rows=20):
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


kv_pairs = [{"label":name, "value":val} for val, name in enumerate(star_names_list)]
vk_pairs = {val:name for val, name in enumerate(star_names_list)}

### Make a dashboard
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = Dash(__name__, external_stylesheets=external_stylesheets)
server = app.server

app.layout = html.Div([
    html.Div(
        children=[
            html.H1("Supernova Arrival Time Delays", style={'text-align': 'center'}),
            dcc.Dropdown(id="candid_selected", options=kv_pairs, multi=False, value=123,
                         placeholder="Select a Star to Explode", style={'width': "90%"}),
            dcc.DatePickerSingle(id='my-date-picker-single', min_date_allowed=date(2022, 4, 14),
                                 max_date_allowed=date(2030, 12, 12), initial_visible_month=date(2022, 6, 14),
                                 date=date(2022, 6, 14)),
            html.Div(id='output-container-date-picker-single')
        ],
        style={'display': 'inline-block', 'width': "30%", 'vertical-align': 'top', 'margin-left': '3vw',
               'margin-top': '3vw'}
    ),

    html.Div(
        children=[
            html.Div([html.H4(children='Delays'),
                      html.Table(id='delay_df', style={"font-weight": "bold", "color": "white"})])
        ],
        style={'display': 'inline-block', 'vertical-align': 'top', 'margin-left': '3vw', 'margin-top': '3vw'}
    ),
    html.Br(),
    html.Div(children=[
        dcc.Graph(id='my_3dscat', figure={})
    ],
        style={'display': 'inline-block', 'vertical-align': 'top', 'margin-left': '3vw', 'margin-top': '3vw'}),
],
    style={'display': 'flex', 'backgroundColor': 'black'})


# ------------------------------------------------------------------------------
# Connect the Plotly graphs with Dash Components
@app.callback(
    [Output(component_id='delay_df', component_property='children'),
     Output(component_id='my_3dscat', component_property='figure'),
     Output(component_id='output-container-date-picker-single', component_property='children')],
    [Input(component_id='my-date-picker-single', component_property='date'),
     Input(component_id='candid_selected', component_property='value')]
)
def update_graph(obstime, candid_selected):
    print(candid_selected)
    star_name = vk_pairs[candid_selected]
    if obstime is not None:
        date_object = date.fromisoformat(obstime)
        date_string = date_object.strftime('%Y-%m-%d')+" 20:00:00.100"
    else:
        date_string = "2022-06-14 20:00:00.100"
    all_delays = get_delays_for_time(obstime_str=date_string)
    selected_df = get_delays_for_star(star_name, delays_dict=all_delays)
    fig = plot_spheres(obstime, candid_name=star_name)
    selected_df['Detectors'] = selected_df.index
    selected_df = selected_df[selected_df.columns.tolist()[::-1]]
    return generate_table(selected_df), fig, date_string


# ------------------------------------------------------------------------------
if __name__ == '__main__':
    app.run_server(debug=True)

