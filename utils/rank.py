def get_overlapped_list(year_distance_threshold=2):
    # dv = {}
    # dl = {}
    # for site in SITES:
    #     vegetation_structure = list_available_urls(VEGETATION_STRUCTURE, site)
    #     lidar = list_available_urls(LIDAR, site)

    #     vegetation_structure = [i.split('/')[-1] for i in vegetation_structure]
    #     lidar = [i.split('/')[-1] for i in lidar]
    #     dv[site] = {}
    #     dl[site] = {}

    #     for v in vegetation_structure:
    #         y, m = v.split('-')
    #         y = int(y)
    #         m = int(m)
    #         if y not in dv[site]:
    #             dv[site][y] = [m]
    #         else:
    #             dv[site][y].append(m)
    #     for v in lidar:
    #         y, m = v.split('-')
    #         y = int(y)
    #         m = int(m)
    #         if y not in dl[site]:
    #             dl[site][y] = [m]
    #         else:
    #             dl[site][y].append(m)
    import json
    with open('/home/toanngo/Documents/GitHub/prisma/preprocessing/dv.json', 'r') as f:
        dv = json.load(f)
    with open('/home/toanngo/Documents/GitHub/prisma/preprocessing/dl.json', 'r') as f:
        dl = json.load(f)
    ranks = {}
    for site in SITES:
        dl_site = dl[site]
        dv_site = dv[site]
        ranks[site] = {}
        for year, dv_months in dv_site.items():
            to_check_dl_years = [str(y) for y in range(int(year) - year_distance_threshold,
                                                       int(year) + year_distance_threshold)]
            scores = []
            infos = []
            for y in to_check_dl_years:
                if y in dl_site:
                    score = 99
                    lidar_date = None
                    dl_months = dl_site[y]
                    min_dv_month = min(dv_months)
                    max_dv_month = max(dv_months)
                    for dl_m in dl_months:
                        if min_dv_month <= dl_m <= max_dv_month:
                            r = 0
                            dl_m_r = dl_m
                        else:
                            r = min(abs(dl_m-min_dv_month),
                                    abs(dl_m-max_dv_month))
                            dl_m_r = dl_m
                        if r < score:
                            score = r
                            lidar_date = f'{dl_m_r}-{y}'
                    cross_year_score = abs(int(y) - int(year))/10.0
                    score += cross_year_score
                    scores.append(score)
                    infos.append(lidar_date)
            if scores:
                min_idx = np.argmin(scores)
                ranks[site][year] = {
                    'score': scores[min_idx],
                    'lidar_date': infos[min_idx],
                    'dv_months': dv_months
                }
            else:
                ranks[site][year] = {}
    return dv, dl, ranks