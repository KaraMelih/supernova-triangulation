import json

base_message = { "False Alarm Prob": "13.96%",
     "_id": "SNEWS_Coincidence_ALERT-UPDATE 2022-09-28T14:38:00.420080",
     "detector_names": ["KM3NeT"],
     "neutrino_times": ["2022-06-14T20:00:00.008268"],
     "p_values": [0.1],
     "p_values average": 0.09285714285714286,
     "sent_time": "2022-09-28T14:38:00.420080",
     "server_tag": "iap-nb-034",
     "sub list number": 0}


def get_json_data(date_string, delay_df):
    message = base_message.copy()
    message['_id'] = f"SNEWS_Coincidence_ALERT-UPDATE {date_string}"
    # date_stripped = ".".join(date_string.strip('.')[:-1])[:-1]
    message["neutrino_times"] = list(delay_df['times'])
    message["sent_time"] = date_string
    message["detector_names"] = list(delay_df.index)
    data_string = json.dumps(message)
    href = f"data:text/json;charset=utf-8,{data_string}"
    return href


base_detector_msg = {
    "detector_name":
        "detector name",
    "neutrino_time":
        "2022-09-28T01:23:45:059496",
    "testing":
        "this is a test"
}

def get_json_per_detector(delay_df):
    master_dict = {}
    for detector, delay in zip(delay_df.index, delay_df['times']):
        message = base_detector_msg.copy()
        message["neutrino_time"] = delay
        message["detector_name"] = detector
        master_dict[detector] = message
    data_string = json.dumps(master_dict)
    href = f"data:text/json;charset=utf-8,{data_string}"
    return href
