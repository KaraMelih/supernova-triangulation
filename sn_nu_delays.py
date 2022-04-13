
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
import numpy as np
import pandas as pd
import plotly
import plotly.tools as tls
import chart_studio.plotly as py
import plotly.graph_objs as go
import plotly.express as px
from plotly.offline import iplot
from dash import Dash, dcc, html, Input, Output

def get_xyz_from_gcrs(gcrs_coor, radius=10, scale_rad=1):
    """ radius in kpc
        returns x,y,z in GCRS coor, in meters
    """
    delta = gcrs_coor.dec.rad
    phi = gcrs_coor.ra.rad
    SNx = -np.cos(phi) * np.cos(delta)
    SNy = -np.sin(phi) * np.cos(delta)
    SNz = -np.sin(delta)
    radius = (radius*u.kpc).to(u.m).value * scale_rad
    return radius * np.array([ SNx, SNy, SNz ])

EarthRadius = 6.3781e6 # in units m
Earthradkpc = EarthRadius*u.m.to(u.kpc)

radius = EarthRadius # in units m
uu, vv = np.mgrid[0:2*np.pi:200j, 0:np.pi:100j]
xE = radius * np.cos(uu)*np.sin(vv)
yE = radius * np.sin(uu)*np.sin(vv)
zE = radius * np.cos(vv)

# delays
df = pd.read_csv('./assets/delays_w-candids.csv')
df.rename(columns={"Unnamed: 0": "Star Name", "Unnamed: 1": "Detector"}, inplace=True)
df.set_index("Star Name", inplace=True)
candid_stars = np.unique(df.index)

# detectors
detectors = pd.read_csv("./assets/detector_locations-xyz.csv")

# candid star positions
radec = np.loadtxt('./assets/star_candidates.txt', usecols=(1, 2), delimiter=',')
ra = radec[:, 0]
dec = radec[:, 1]
timestring = "2022-06-14 20:00:00.100"
ObsTime = Time(timestring, format='iso', scale='utc')
skycoors = SkyCoord(ra, dec, frame="fk5", obstime=ObsTime, unit='deg')
candid_xyz = get_xyz_from_gcrs(skycoors.gcrs,
                               radius=Earthradkpc * np.random.uniform(low=1.5, high=12, size=len(skycoors)))

candid_names = np.loadtxt('./assets/star_candidates.txt', usecols=0, delimiter=',', dtype=str)
_names = []
for name in candid_names:
    stripped_name = name.strip()
    if stripped_name not in _names:
        _names.append(stripped_name)
    else:
        _names.append(stripped_name + "-2")

df_candids = pd.DataFrame(candid_xyz.T, columns=('x', 'y', 'z'))
df_candids['RA'] = ra
df_candids['DEC'] = dec
df_candids['names'] = _names
df_candids.set_index('names', inplace=True)

candid_text = [name + f"<br>RA:{ra:.2f}" + f"<br>DEC:{dec:.2f}" for name, ra, dec in zip(df_candids.index, df_candids['RA'], df_candids['DEC'])]
detector_text = [name + f"<br>lon: {lon:.2f}deg<br>lat: {lat:.2f}deg<br>Height{h:.2f}m" for name,lon,lat,h in zip(detectors["Unnamed: 0"],
                                                                                                                  detectors["Longitude (deg)"],
                                                                                                                  detectors["Latitude (deg)"],
                                                                                                                  detectors["Height (m)"])]

lines = []
line_marker = dict(color='green', width=3) # '#0066FF'
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

td_data =[go.Scatter3d(x=df_candids['x'], y=df_candids['y'], z=df_candids['z'], mode ='markers', hoverinfo='text', text=candid_text, marker=dict(size=5))]
td_data2 =[go.Scatter3d(x=df_candids['x'], y=df_candids['y'], z=df_candids['z'], mode ='markers', hoverinfo='none', marker=dict(size=10, color='white', opacity=0.5))]
det_data =[go.Scatter3d(x=detectors['x (m)'], y=detectors['y (m)'], z=detectors['z (m)'], marker=dict(size=7, symbol='x'), text=detector_text, hoverinfo='text', mode ='markers')]

def plot_spheres(candid_name=None):
    data = lines+td_data2+td_data+det_data
    if candid_name is not None:
        candid = df_candids.loc[candid_name]
        print(candid)
        modified_part =[go.Scatter3d(x=[candid['x']], y=[candid['y']], z=[candid['z']],
                                     mode='markers', hoverinfo='text', marker=dict(size=20, color='yellow', opacity=0.7),
                                     text=f"{candid_name} <br> EXPLODED!")]
        data = data + modified_part
    fig = go.Figure(data=data, layout=layout)
#     iplot(fig, filename='wireframe_plot')
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


kv_pairs = [{"label":name, "value":val} for val, name in enumerate(candid_stars)]
vk_pairs = {val:name for val, name in enumerate(candid_stars)}

### Make a dashboard
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = Dash(__name__, external_stylesheets=external_stylesheets)
server = app.server

app.layout = html.Div([
    html.Div(children=[
                html.H1("Supernova Arrival Time Delays", style={'text-align': 'center'}),
                dcc.Dropdown(id="candid_selected",
                             options=kv_pairs,
                             multi=False,
                             value=123,
                             placeholder="Select a Star to Explode",
                             style={'width': "90%"}
                             )],
            style={'display': 'inline-block','width': "30%", 'vertical-align': 'top', 'margin-left': '3vw', 'margin-top': '3vw'}
    ),

    html.Div(children=[
                html.Div([html.H4(children='Delays'), html.Table(id='delay_df',
                                                                 style={"font-weight": "bold", "color":"white"})],
                         )],
        style={'display': 'inline-block', 'vertical-align': 'top', 'margin-left': '3vw', 'margin-top': '3vw'}
    ),
    html.Br(),
    html.Div(children=[
                    dcc.Graph(id='my_3dscat', figure={})
    ],
        style={'display': 'inline-block', 'vertical-align': 'top', 'margin-left': '3vw', 'margin-top': '3vw'}),
],
        style={'display': 'flex', 'backgroundColor':'black'})


# ------------------------------------------------------------------------------
# Connect the Plotly graphs with Dash Components
@app.callback(
    [Output(component_id='delay_df', component_property='children'),
     Output(component_id='my_3dscat', component_property='figure')],
    [Input(component_id='candid_selected', component_property='value')]
)
def update_graph(candid_selected):
    print(candid_selected)
    print(type(candid_selected))
    name = vk_pairs[candid_selected]
    print(f">>> |{name}|")
    # container = f"><><{name}><><" #df.loc[name] #"The year chosen by user was: {}".format(option_slctd)
    fig = plot_spheres(name)
    return  generate_table(df.loc[name]), fig


# ------------------------------------------------------------------------------
if __name__ == '__main__':
    app.run_server(debug=True)



