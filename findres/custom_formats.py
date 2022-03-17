from obspy import UTCDateTime

def _rise_phases(origin_time, station, phase_file):
    p_pick, s_pick = None, None
    with phase_file.open('r') as file:
        start_section = False
        event_origin = None
        while line := file.readline():
            if line.isspace():
                if start_section:
                    break
                else:
                    continue
            fields = line.split()
            if len(fields) == 12:
                event_origin = fields[2] + fields[3] + fields[4] + fields[7] + fields[8] + fields[9]
            if event_origin and event_origin.startswith(origin_time.strftime('%Y%m%d%H%M%S')):
                start_section = True
                pick_datetime = fields[6] + fields[7] + fields[8]
                pick_station = fields[0]
                pick_type = fields[4]
                if pick_station == station:
                    if pick_type == 'P':
                        p_pick = UTCDateTime.strptime(pick_datetime, '%Y%m%d%H%M%S.%f')
                    elif pick_type == 'S':
                        s_pick = UTCDateTime.strptime(pick_datetime, '%Y%m%d%H%M%S.%f')
    
    return p_pick, s_pick

def _rise_errors(origin_time, phase_file):
    return None 

formats_registry = {'rise_custom': {'phases': _rise_phases, 'errors': _rise_errors}}