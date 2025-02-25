def get_addresses(groups, i, s, d, lim):
    """
    Updates the dictionaries i, s, d and the counter lim based on the groups.
    
    Parameters:
        groups (list of lists): Each element is a list of strings.
        i (dict): Nested dictionary that will store the numeric addresses.
        s (dict): Dictionary mapping each group element (string) to a list of addresses.
        d (dict): Dictionary mapping an integer (lim) to the corresponding address as a list of strings.
        lim (int): Starting counter for addresses.
    
    Returns:
        tuple: (i, s, d, lim) updated with the new addresses.
    """
    # Ensure that every element in groups exists in s (initialize to empty list if not)
    fnames = list(s.keys()) if s else []
    for gp in groups:
        for elem in gp:
            if elem not in fnames:
                s[elem] = []

    n = len(groups)
    
    if n == 1:
        gp1 = groups[0]
        for elem in gp1:
            lim += 1
            # Assign directly since there's only one level.
            i[elem] = lim
            s[elem].append(lim)
            d[lim] = [elem]
            
    elif n == 2:
        gp1, gp2 = groups[0], groups[1]
        for elem1 in gp1:
            # Ensure a nested dictionary exists for elem1.
            if elem1 not in i or not isinstance(i[elem1], dict):
                i[elem1] = {}
            for elem2 in gp2:
                lim += 1
                i[elem1][elem2] = lim
                s[elem1].append(lim)
                s[elem2].append(lim)
                d[lim] = [elem1, elem2]
                
    elif n == 3:
        gp1, gp2, gp3 = groups[0], groups[1], groups[2]
        for elem1 in gp1:
            if elem1 not in i or not isinstance(i[elem1], dict):
                i[elem1] = {}
            for elem2 in gp2:
                if elem2 not in i[elem1] or not isinstance(i[elem1][elem2], dict):
                    i[elem1][elem2] = {}
                for elem3 in gp3:
                    lim += 1
                    i[elem1][elem2][elem3] = lim
                    s[elem1].append(lim)
                    s[elem2].append(lim)
                    s[elem3].append(lim)
                    d[lim] = [elem1, elem2, elem3]
                    
    elif n == 4:
        gp1, gp2, gp3, gp4 = groups[0], groups[1], groups[2], groups[3]
        for elem1 in gp1:
            if elem1 not in i or not isinstance(i[elem1], dict):
                i[elem1] = {}
            for elem2 in gp2:
                if elem2 not in i[elem1] or not isinstance(i[elem1][elem2], dict):
                    i[elem1][elem2] = {}
                for elem3 in gp3:
                    if elem3 not in i[elem1][elem2] or not isinstance(i[elem1][elem2][elem3], dict):
                        i[elem1][elem2][elem3] = {}
                    for elem4 in gp4:
                        lim += 1
                        i[elem1][elem2][elem3][elem4] = lim
                        s[elem1].append(lim)
                        s[elem2].append(lim)
                        s[elem3].append(lim)
                        s[elem4].append(lim)
                        d[lim] = [elem1, elem2, elem3, elem4]
                        
    elif n == 5:
        gp1, gp2, gp3, gp4, gp5 = groups[0], groups[1], groups[2], groups[3], groups[4]
        for elem1 in gp1:
            if elem1 not in i or not isinstance(i[elem1], dict):
                i[elem1] = {}
            for elem2 in gp2:
                if elem2 not in i[elem1] or not isinstance(i[elem1][elem2], dict):
                    i[elem1][elem2] = {}
                for elem3 in gp3:
                    if elem3 not in i[elem1][elem2] or not isinstance(i[elem1][elem2][elem3], dict):
                        i[elem1][elem2][elem3] = {}
                    for elem4 in gp4:
                        if elem4 not in i[elem1][elem2][elem3] or not isinstance(i[elem1][elem2][elem3][elem4], dict):
                            i[elem1][elem2][elem3][elem4] = {}
                        for elem5 in gp5:
                            lim += 1
                            i[elem1][elem2][elem3][elem4][elem5] = lim
                            s[elem1].append(lim)
                            s[elem2].append(lim)
                            s[elem3].append(lim)
                            s[elem4].append(lim)
                            s[elem5].append(lim)
                            d[lim] = [elem1, elem2, elem3, elem4, elem5]

    # Record the total number of states.
    i["nstates"] = lim

    return i, s, d, lim
