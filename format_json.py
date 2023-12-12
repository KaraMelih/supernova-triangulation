import json
# import snews_pt

base_alert_message = {"False Alarm Prob": "0.00%",
                      "_id": "SNEWS_Coincidence_ALERT datestring",
                      "alert_type": "COINC_MSG",
                      "detector_names": ["IceCube", "NOvA"],
                      "neutrino_times": ["2030-01-01T12:30:19.678999", "2030-01-01T12:30:21.678999" ],
                      "p_values": [0.98, 0.98],
                      "p_values average": 0.98,
                      "sent_time": "2023-12-12T12:42:58.225935",
                      "server_tag": "iap-nb-034",
                      "sub list number": 0}

def get_json_data(date_string, delay_df):
    message = base_alert_message.copy()
    message['_id'] = f"SNEWS_Coincidence_ALERT {date_string}"
    message["neutrino_times"] = list(delay_df['times'])
    message["sent_time"] = date_string
    message["detector_names"] = list(delay_df.index)
    data_string = json.dumps(message)
    href = f"data:text/json;charset=utf-8,{data_string}"
    return href


base_detector_msg = {"detector_name": "detector name",
                     "neutrino_time": "2022-09-28T01:23:45:059496",
                     "is_test": True}

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
