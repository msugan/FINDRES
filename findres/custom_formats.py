from obspy import UTCDateTime

def _phs_phases(origin_time, station, phase_file):
    p_pick, s_pick = None, None
    with phase_file.open('r') as file:
        start_section = False
        while line := file.readline():
            fields = line.split()
            entry_station = fields[0]
            entry_pick = fields[4]
            entry_datetime = fields[6] + fields[7] + fields[8]
            if entry_datetime.startswith(origin_time.strftime('%Y%m%d%H%M%S')):
                start_section = True
                if entry_station == station:
                    if entry_pick == 'P':
                        p_pick = UTCDateTime.strptime(entry_datetime, '%Y%m%d%H%M%S%f')
                    elif entry_pick == 'S':
                        s_pick = UTCDateTime.strptime(entry_datetime, '%Y%m%d%H%M%S%f')
            if start_section and line.isspace():
                break
    
    return p_pick, s_pick

def _phs_errors(origin_time, phase_file):
    return None 

formats_registry = {'phs': {'phases': _phs_phases, 'errors': _phs_errors}}