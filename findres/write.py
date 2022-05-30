from itertools import combinations


def dump_hypodd(families, catalogue, errors, repeaters, parameters, output_dir):
    for i, family in enumerate(families):
        with (output_dir / f"{i}_RES_event.sel").open("w") as file:
            for n in family:
                e_t, e_h, e_v = map(lambda s: err if (err := errors[catalogue.loc[n, "name"]].get(s)) else parameters[
                    'hypodd_default_' + s], ['time_uncertainty', 'horizontal_uncertainty', 'vertical_uncertainty'])
                file.write(f"{catalogue.loc[n, 'date'].strftime('%Y%m%d  %H%M%S%f')[:-4]}   "
                           f"{catalogue.loc[n, 'latitude']:.4f}     {catalogue.loc[n, 'longitude']:.4f}    "
                           f"{catalogue.loc[n, 'depth']:.3f}   {catalogue.loc[n, 'magnitude']:.1f}    "
                           f"{e_t:.2f}    {e_h:.2f}   {e_v:.2f}        "
                           f"{catalogue.loc[n, 'name']}\n")
        if len(family) > 2:
            with (output_dir / f"{i}_RES_dt.cc").open("w") as file:
                for (t1, t2) in combinations(family, 2):
                    if repeaters[(t1, t2)]:
                        file.write(f"#    {catalogue.loc[t1, 'name']}    {catalogue.loc[t2, 'name']}     0.0\n")
                        for (network, station) in sorted(repeaters[(t1, t2)]):
                            cc, delta_sp = repeaters[(t1, t2)][(network, station)]
                            delta_v = parameters['Vp'] - parameters['Vs']
                            ttp = -parameters['Vs'] * delta_sp / delta_v
                            file.write(f"{station}     {ttp: 10.9f}    {cc:.2f}    P\n")
                            tts = -parameters['Vp'] * delta_sp / delta_v
                            file.write(f"{station}     {tts: 10.9f}    {cc:.2f}    S\n")
