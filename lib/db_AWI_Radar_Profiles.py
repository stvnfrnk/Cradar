
def get_AWI_radar_profile_metadata():

    import numpy as np
    import pandas as pd
    import os, glob

    ###################
    # All ASIRAS Season

    dict_ASIRAS    = {"Campaign Folder"  : ["CryoVEx_2004B",
                                            "CryoVEx_2004C",
                                            "CryoVEx_2005A",
                                            "CryoVEx_2006A",
                                            "CryoVEx_2006B",
                                            "CryoVEx_2007A",
                                            "CryoVEx_2007B",
                                            "CryoVEx_2007C",
                                            "CryoVEx_2007D",
                                            "CryoVEx_2008A",
                                            "CryoVEx_2008B",
                                            "CryoVEx_2008C",
                                            "CryoVEx_2009A",
                                            "CryoVEx_2011A",
                                            "CryoVEx_2011B",
                                            "CryoVEx_2012A",
                                            "CryoVEx_2014",
                                            "CryoVEx_2016A",
                                            "CryoVEx_2016B",
                                            "CryoVEx_2017",
                                            "CryoVEx_ANT2014",
                                            "CryoVEx_ANT2016",
                                            "DTU_2016A",
                                            "DTU_2016B",
                                            "DTU_2017",
                                            "DTU_2019",
                                            "DTU_ANT2018"],

                    "Season Name"       : ["ARK 2004",
                                            "ARK 2004",
                                            "ARK 2005",
                                            "ARK 2006",
                                            "ARK 2006",
                                            "ARK 2007",
                                            "ARK 2007",
                                            "ARK 2007",
                                            "ANT 2007/08",
                                            "ARK 2008",
                                            "ARK 2008",
                                            "ANT 2008/09",
                                            "ARK 2009",
                                            "ARK 2011",
                                            "ANT 2011/12",
                                            "ARK 2012",
                                            "ARK 2014",
                                            "ARK 2016",
                                            "ARK 2016",
                                            "ARK 2017",
                                            "ANT 2014/15",
                                            "ANT 2015/16",
                                            "ARK 2016",
                                            "ARK 2016",
                                            "ARK 2017",
                                            "ARK 2019",
                                            "ANT 2017/18"],

                    "Radar System"     : "ASIRAS",
                    "Radar Specs"    : "13-14.5 GHz bandwidth",

                    "Platform"         :  ["tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd"],
                                        
                    "Campaign Name"    :  ["CryoVEx_2004B",
                                            "CryoVEx_2004C",
                                            "CryoVEx_2005A",
                                            "CryoVEx_2006A",
                                            "CryoVEx_2006B",
                                            "CryoVEx_2007A",
                                            "CryoVEx_2007B",
                                            "CryoVEx_2007C",
                                            "CryoVEx_2007D",
                                            "CryoVEx_2008A",
                                            "CryoVEx_2008B",
                                            "CryoVEx_2008C",
                                            "CryoVEx_2009A",
                                            "CryoVEx_2011A",
                                            "CryoVEx_2011B",
                                            "CryoVEx_2012A",
                                            "CryoVEx_2014",
                                            "CryoVEx_2016A",
                                            "CryoVEx_2016B",
                                            "CryoVEx_2017",
                                            "CryoVEx_ANT2014",
                                            "CryoVEx_ANT2016",
                                            "DTU_2016A",
                                            "DTU_2016B",
                                            "DTU_2017",
                                            "DTU_2019",
                                            "DTU_ANT2018"],

                    "Campaign Name Long":  ["CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "CryoSat-2 Validation Experiment",
                                            "DTU_2016A",
                                            "DTU_2016B",
                                            "DTU_2017",
                                            "DTU_2019",
                                            "DTU_ANT2018"],

                    "Campaign PI"       : np.repeat("None", 27),
                    "Grant ID"          : np.repeat("None", 27),
                    "Pangaea Identifyer" : np.repeat("None", 27),
                    "Radar Data DOI"     : np.repeat("None", 27)}

    ## ================================================================================================= ##
    ## ================================================================================================= ##



    ###################
    # All SNOW Season

    dict_SNOW     = {"Campaign Folder"  : ["antr2014",
                                        "antr2015",
                                        "arkr2014",
                                        "arkr2015"],

                    "Season Name"      : ["ANT 2013/14",
                                        "ANT 2014/15",
                                        "ARK 2014",
                                        "ARK 2015"],

                    "Radar System"     : "SNOW",
                    "Radar Specs"    : "4-8 GHz bandwidth",

                    "Platform"         : ["Polar 6",
                                        "Polar 6",
                                        "Polar 6",
                                        "Polar 6"],
                                        
                    "Campaign Name"    : ["CofiStructure",
                                        "tbd", 
                                        "tbd", 
                                        "tbd"],

                    "Campaign Name Long": ["tbd",
                                        "tbd", 
                                        "tbd", 
                                        "tbd"],

                    "Campaign PI"      : ["Daniel Steinhage (AWI)",
                                        "Daniel Steinhage (AWI)",
                                        "Daniel Steinhage (AWI)",
                                        "Daniel Steinhage (AWI)"],

                    "Grant ID"        : ["None",
                                        "None", 
                                        "None", 
                                        "None"],

                    "Pangaea Identifyer" : ["tbd",
                                            "tbd",
                                            "tbd",
                                            "tbd"],

                    "Radar Data DOI"     : ["https://doi.org/10.1594/PANGAEA.983487",
                                            "None",
                                            "None",
                                            "None"]}

    ## ================================================================================================= ##
    ## ================================================================================================= ##


    ###################
    # All ACCU Season

    dict_ACCU     = {"Campaign Folder"  : ["antr2011",
                                        "antr2012",
                                        "antr2013",
                                        "antr2014",
                                        "antr2015",
                                        "antr2017",
                                        "arkr2010",
                                        "arkr2012",
                                        "arkr2015",
                                        "arkr2016",
                                        "arkr2018"],

                    "Season Name"      : ["ANT 2010/11",
                                        "ANT 2011/12",
                                        "ANT 2012/13",
                                        "ANT 2013/14",
                                        "ANT 2014/15",
                                        "ANT 2016/17",
                                        "ARK 2010",
                                        "ARK 2012",
                                        "ARK 2015",
                                        "ARK 2016",
                                        "ARK 2018"],

                    "Radar System"     : "ACCU",
                    
                    "Radar Specs"    : ["400–800 MHz bandwidth",
                                        "400–800 MHz bandwidth",
                                        "400–800 MHz bandwidth",
                                        "400–800 MHz bandwidth",
                                        "400–800 MHz bandwidth",
                                        "400–800 MHz bandwidth",
                                        "500–700 MHz bandwidth",
                                        "400–800 MHz bandwidth",
                                        "400–800 MHz bandwidth",
                                        "400–800 MHz bandwidth",
                                        "400–800 MHz bandwidth"],

                    "Platform"         : ["Polar 6",
                                        "Polar 6",
                                        "Polar 6",
                                        "Polar 6",
                                        "Polar 6",
                                        "Polar 6",
                                        "Polar 5",
                                        "Polar 5",
                                        "Polar 6",
                                        "Polar 6",
                                        "Polar 6"],
                                        
                    "Campaign Name"    : ["WEGAS",
                                        "GEA II, WEGAS", 
                                        "WEGAS, MaBaJu 1, GEA III, RECISL", 
                                        "CoFi Structure, MaBaJu 2, RECISL, Sör-Mond, VELMAP", 
                                        "GEA IV", 
                                        "GEA-Vb-FMA, OIR Dome F",
                                        "NEEM & EGIG", 
                                        "Structure NGT 2012",
                                        "MABANG 1",
                                        "tbd",
                                        "tbd"],

                    "Campaign Name Long": ["West-East Gondwana Amalgamation and its Separation",
                                        "GEA II: Geodynamic Evolution of East Antarctica (Part 2), WEGAS: West-East Gondwana Amalgamation and its Separation", 
                                        "WEGAS: West-East Gondwana Amalgamation and its Separation, MaBaJu 1: Mass balance of Jutulstraumen (Part 1), GEA III: Geodynamic Evolution of East Antarctica (Part 3), RECISL: Recovery Glacier and Lakes, ice thickness, bedrock topography and basal properties", 
                                        "CoFi Structure, MaBaJu 2, RECISL, Sör-Mond, VELMAP", 
                                        "Geodynamic Evolution of East Antarctica (Part 4)", 
                                        "GEA-Vb-FMA: Geodynamic Evolution of East Antarctica, the Forster Magnetic Anomaly, OIR Dome F: Oldest Ice Reconnaissance Dome Fuji",
                                        "North EEMian ice core project and International Glaciological Expedition to Greenland (EGIG) line", 
                                        "Structure of the Greenland ice sheet along the 77° North Greenland Traverse 2012",
                                        "Mass Balance of Glaciers in North Greenland (Part 1)",
                                        "tbd",
                                        "tbd"],

                    "Campaign PI"      : ["Daniel Steinhage (AWI)",
                                        "Daniel Steinhage (AWI)",
                                        "Daniel Steinhage (AWI)",
                                        "Daniel Steinhage (AWI)",
                                        "Daniel Steinhage (AWI)",
                                        "Daniel Steinhage (AWI)",
                                        "Heinrich Miller (AWI)",
                                        "Heinrich Miller (AWI)",
                                        "Angelika Humbert (AWI)",
                                        "Daniel Steinhage (AWI)",
                                        "Daniel Steinhage (AWI)"],

                    "Grant ID"        :  ["AWI_PA_02020",
                                        "AWI_PA_02030, AWI_PA_02031", 
                                        "AWI_PA_02039, AWI_PA_02040, AWI_PA_02041, AWI_PA_02043", 
                                        "AWI_PA_02037, AWI_PA_02038, AWI_PA_02056, AWI_PA_02058", 
                                        "AWI_PA_02068", 
                                        "AWI_PA_02083, AWI_PA_02084",
                                        "AWI_PA_02034", 
                                        "AWI_PA_02017",
                                        "AWI_PA_02072",
                                        "None",
                                        "None"],

                    "Pangaea Identifyer" : ["WEGAS_2010/11",
                                            "WEGAS_2011/12", 
                                            "WEGAS_2012/13", 
                                            "ANT_2013/14",
                                            "ANT_2014/15", 
                                            "P6-204_ANT_OIR_2016_2017",
                                            "NEEM_2010",
                                            "NGT_2012",
                                            "P6-196_MABANG_2015",
                                            "tbd",
                                            "tbd"],

                    "Radar Data DOI"     : ["https://doi.org/10.1594/PANGAEA.983392",  # antr2011
                                            "https://doi.org/10.1594/PANGAEA.983396",  # antr2012
                                            "https://doi.org/10.1594/PANGAEA.983435",  # antr2013
                                            "https://doi.org/10.1594/PANGAEA.983439",  # antr2014
                                            "https://doi.org/10.1594/PANGAEA.983468",  # antr2015
                                            "https://doi.org/10.1594/PANGAEA.983479",  # antr2017
                                            "https://doi.org/10.1594/PANGAEA.982840",  # arkr2010
                                            "https://doi.org/10.1594/PANGAEA.982845",  # arkr2012
                                            "https://doi.org/10.1594/PANGAEA.982772",  # arkr2015
                                            "None",  # arkr2016 (sea ice)
                                            "None"]  # arkr2018 (sea ice)
                                        }

    ## ================================================================================================= ##
    ## ================================================================================================= ##


    ###################
    # All EMR Season

    # ANT
    dict_EMR_ANT  = {"Campaign Folder"  : ["antr1995",
                                        "antr1996",
                                        "antr1997",
                                        "antr1998",
                                        "antr1999",
                                        "antr2001",
                                        "antr2002",
                                        "antr2003",
                                        "antr2004",
                                        "antr2005",
                                        "antr2006",
                                        "antr2008",
                                        "antr2009",
                                        "antr2011",
                                        "antr2012",
                                        "antr2013",
                                        "antr2014",
                                        "antr2015",
                                        "antr2016",
                                        "antr2017",
                                        "antr2018",
                                        "antr2023"],

                    "Season Name"      : ["ANT 1994/95",
                                        "ANT 1995/96",
                                        "ANT 1996/97",
                                        "ANT 1997/98",
                                        "ANT 1998/99",
                                        "ANT 2000/01",
                                        "ANT 2001/02",
                                        "ANT 2002/03",
                                        "ANT 2003/04",
                                        "ANT 2004/05",
                                        "ANT 2005/06",
                                        "ANT 2007/08",
                                        "ANT 2008/09",
                                        "ANT 2010/11",
                                        "ANT 2011/12",
                                        "ANT 2012/13",
                                        "ANT 2013/14",
                                        "ANT 2014/15",
                                        "ANT 2015/16",
                                        "ANT 2016/17",
                                        "ANT 2017/18",
                                        "ANT 2022/23"],

                    "Radar System"     : "EMR",
                    "Radar Specs"    : "from ID",

                    "Platform"         : ["Polar 2",
                                        "Polar 2",
                                        "Polar 2",
                                        "Polar 2",
                                        "Polar 2",
                                        "Polar 2",
                                        "Polar 2",
                                        "Polar 2",
                                        "Polar 2",
                                        "Polar 2",
                                        "Polar 2",
                                        "Polar 5",
                                        "Polar 5",
                                        "Polar 5",
                                        "Polar 6",
                                        "Polar 6",
                                        "Polar 6",
                                        "Polar 6",
                                        "Polar 5",
                                        "Polar 6",
                                        "Polar 6",
                                        "Polar 5"],

                    "Campaign Name"    : ["MAGRAD",                            # antr1995
                                        "EPICA I",                           # antr1996
                                        "EPICA II",                          # antr1997
                                        "EPICA III",                         # antr1998
                                        "EPICA IV",                          # antr1999
                                        "EPICA VI, SEAL, VISA",              # antr2001
                                        "EPICA VII",                         # antr2002
                                        "EPICA VIII, VISA II, SEAL",         # antr2003
                                        "EPICA IX, VISA III, SEAL",          # antr2004
                                        "EPICA X, VISA IV",                  # antr2005
                                        "ANTSYO",                            # antr2006
                                        "DoCo, VISA",                        # antr2008
                                        "RBI, WEGAS",                        # antr2009
                                        "GEA I, WEGAS, EUFAR",               # antr2011
                                        "GEA II, WEGAS",                     # antr2012
                                        "WEGAS, MaBaJu 1, GEA III, RECISL",  # antr2013
                                        "WEGAS, SWIFT, STRUCGLAC, MaBaJu 2, RECISL, Sör-Mond, VELMAP", # antr2014
                                        "GEA IV",                            # antr2015
                                        "GEA-Va-FMA",                        # antr2016
                                        "GEA-Vb-FMA, OIR Dome F, MaBaJu 2",  # antr2017
                                        "ANIRES",                            # antr2018
                                        "RIISERBATHY"],                      # antr2023

                    "Campaign Name Long": ["MAGRAD",                           # antr1995
                                        "EPICA pre-site survey (Part 1)",                           # antr1996
                                        "EPICA pre-site survey (Part 2)",                          # antr1997
                                        "EPICA pre-site survey (Part 3)",                         # antr1998
                                        "EPICA pre-site survey (Part 4)",                          # antr1999
                                        "EPICA VI: EPICA pre-site survey (Part 6), SEAL, VISA: Validierung, Verdichtung und Interpretation von Satellitendaten zur Bestimmung (Part 1)",
                                        "EPICA VII: EPICA pre-site survey (Part 7)",
                                        "EPICA VIII: EPICA pre-site survey (Part 8), VISA II: Validierung, Verdichtung und Interpretation von Satellitendaten zur Bestimmung (Part 2), SEAL",
                                        "EPICA IX: EPICA pre-site survey (Part 9), VISA III: Validierung, Verdichtung und Interpretation von Satellitendaten zur Bestimmung (Part 3), SEAL",
                                        "EPICA X: EPICA pre-site survey (Part 10), VISA IV: Validierung, Verdichtung und Interpretation von Satellitendaten zur Bestimmung (Part 4)",
                                        "Antarctic Flight Missions at Syowa Region",
                                        "DoCo: Dome Connection East Antarctica, VISA V: Validierung, Verdichtung und Interpretation von Satellitendaten zur Bestimmung (Part 5)",
                                        "RBI: Reconnaissance Berkner Island, WEGAS: West-East Gondwana Amalgamation and its Separation",
                                        "GEA I: Geodynamic Evolution of East Antarctica (Part 1), WEGAS: West-East Gondwana Amalgamation and its Separation, EUFAR: European Facilities for Airborne Research",
                                        "GEA II: Geodynamic Evolution of East Antarctica (Part 2), WEGAS: West-East Gondwana Amalgamation and its Separation",
                                        "WEGAS: West-East Gondwana Amalgamation and its Separation, MaBaJu 1: Mass balance of Jutulstraumen (Part 1), GEA III: Geodynamic Evolution of East Antarctica (Part 3), RECISL: Recovery Glacier and Lakes, ice thickness, bedrock topography and basal properties",
                                        "WEGAS, SWIFT, STRUCGLAC, MaBaJu 2, RECISL, Sör-Mond, VELMAP", # antr2014
                                        "Geodynamic Evolution of East Antarctica (Part 4)",
                                        "Geodynamic Evolution of East Antarctica, the Forster Magnetic Anomaly",
                                        "GEA-Vb-FMA: Geodynamic Evolution of East Antarctica, the Forster Magnetic Anomaly, OIR Dome F: Oldest Ice Reconnaissance Dome Fuji, MaBaJu 2: Mass balance of Jutulstraumen – Part 2",
                                        "Antarctic anisotropy from radio-echo sounding",
                                        "Riiser-Larsen Ice Shelf Geophysical Survey"],

                    "Campaign PI"      : ["Heinrich Miller (AWI)", # antr1995
                                        "Heinrich Miller (AWI)", # antr1996
                                        "Heinrich Miller (AWI)", # antr1997
                                        "Heinrich Miller (AWI)", # antr1998
                                        "Heinrich Miller (AWI)", # antr1999
                                        "Heinrich Miller (AWI)", # antr2001
                                        "Heinrich Miller (AWI)", # antr2002
                                        "Heinrich Miller (AWI)", # antr2003
                                        "Heinrich Miller (AWI)", # antr2004
                                        "Heinrich Miller (AWI)", # antr2005
                                        "Uwe Nixdorf (AWI)",     # antr2006
                                        "Daniel Steinhage (AWI)", # antr2008
                                        "RBI: Olaf Eisen (AWI), WEGAS: Wilfried Jokat (AWI)",       # antr2009
                                        "GEA I & WEGAS: Wilfried Jokat (AWI), EUFAR: tbd", # antr2011
                                        "Wilfried Jokat (AWI)",  # antr2012
                                        "WEGAS & GEA III: Wilfried Jokat (AWI), MaBaJu1 & RECISL: Angelika Humbert (AWI), ", # antr2013
                                        "WEGAS: Wilfried Jokat (AWI), SWIFT: Thomas Kleiner (AWI), STRUCGLAC: Daniel Steinhage (AWI), MaBaJu 2: Angelika Humbert (AWI), RECISL: Angelika Humbert (AWI), Sör-Mond: Olaf Eisen (AWI), VELMAP: Matthias Braun (Unversity of Erlangen)",  # antr2014
                                        "Wilfried Jokat (AWI)", # antr2015
                                        "Graeme Eagles (AWI)", # antr2016
                                        "GEA-Vb-FMA: Graeme Eagles (AWI), OIR Dome F: Olaf Eisen (AWI), MaBaJu 2: Angelika Humbert (AWI)", # antr2017
                                        "Daniel Steinhage (AWI)", #antr2018
                                        "Graeme Eagles (AWI)" #antr2023
                                        ],

                    "Grant ID"          : ["None",    # antr1995
                                        "None",    # antr1996
                                        "None",    # antr1997
                                        "None",    # antr1998
                                        "None",    # antr1999
                                        "None",
                                        "None",
                                        "None",
                                        "None",
                                        "None",
                                        "None",
                                        "DoCo: AWI_PA_02004, WEGAS: AWI_PA_02003",
                                        "RBI: AWI_PA_02009, WEGAS: AWI_PA_02010",
                                        "GEA I: AWI_PA_02021, WEGAS: AWI_PA_02020, EUFAR: tbd",
                                        "GEA II: AWI_PA_02031, WEGAS: AWI_PA_02030",
                                        "GEA III: AWI_PA_02039, WEGAS: AWI_PA_02040, MaBaJu 1: AWI_PA_02043, RECISL: AWI_PA_02041",
                                        "WEGAS: AWI_PA_02051, SWIFT: AWI_PA_02053, STRUCGLAC: AWI_PA_02054, MaBaJu 2: AWI_PA_02055, RECISL: AWI_PA_02056, Sör-Mond: AWI_PA_02058, VELMAP: AWI_PA_02038",
                                        "AWI_PA_02068",
                                        "AWI_PA_02077",
                                        "GEA-Vb-FMA: AWI_PA_02083, OIR Dome F: AWI_PA_02084, MaBaJu 2: AWI_PA_02085",
                                        "AWI_PA_02092",
                                        "AWI_PA_02138"],

                    "Pangaea Identifyer": ["ANT_1994/95",  # antr1995
                                        "EPICA_I",      # antr1996
                                        "EPICA_II",     # antr1997
                                        "EPICA_II",     # antr1998
                                        "EPICA_IV",     # antr1999
                                        "EPICA_VI",     # antr2001
                                        "EPICA_VII",    # antr2002
                                        "EPICA_VIII",   # antr2003
                                        "EPICA_IX",     # antr2004
                                        "EPICA_X",      # antr2005
                                        "ANT_2005/06",  # antr2006
                                        "DoCo_2007/08, VISA_2007, ANT_2007/08", # antr2008
                                        "RBI, WEGAS",    # antr2009
                                        "WEGAS_2010/11", # antr2011
                                        "WEGAS_2011/12, GEA_2012",     # antr2012
                                        "WEGAS_2012/13, MaBaJu_2013",  # antr2013
                                        "ANT_2013/14, GEA-IV_P6_2013", # antr2014
                                        "ANT_2014/15, GEA-IV_P6_2014", # antr2015
                                        "P5-200_ANT_2015_2016",        # antr2016
                                        "P6-204_ANT_GEA_2016_2017, P6-204_ANT_OIR_2016_2017, P6-204_ANT_UWB_2016_2017", # antr2017
                                        "P6-209_ANIRES_2017_2018",     # antr2018
                                        "P5-236_RIISERBATHY_22_23"],   # antr2023

                    "Radar Data DOI"    : ["None", #antr1995
                                        "None", #antr1996
                                        "https://doi.org/10.1594/PANGAEA.974229", #antr1997
                                        "https://doi.org/10.1594/PANGAEA.974295", #antr1998
                                        "https://doi.org/10.1594/PANGAEA.974402", #antr1999
                                        "https://doi.org/10.1594/PANGAEA.974431", #antr2001
                                        "https://doi.org/10.1594/PANGAEA.974432", #antr2002
                                        "https://doi.org/10.1594/PANGAEA.982598", #antr2003
                                        "None", #antr2004
                                        "https://doi.org/10.1594/PANGAEA.983337", #antr2005
                                        "https://doi.org/10.1594/PANGAEA.974441", #antr2006
                                        "https://doi.org/10.1594/PANGAEA.974448", #antr2008
                                        "https://doi.org/10.1594/PANGAEA.974453", #antr2009
                                        "None", #antr2011
                                        "None", #antr2012
                                        "None", #antr2013
                                        "None", #antr2014
                                        "None", #antr2015
                                        "None", #antr2016
                                        "https://doi.org/10.1594/PANGAEA.982745", #antr2017
                                        "None", #antr2018
                                        "https://doi.pangaea.de/10.1594/PANGAEA.971622"]} # antr2023

    # ARK
    dict_EMR_ARK  = {"Campaign Folder"  : ["arkr1995",
                                        "arkr1996",
                                        "arkr1997",
                                        "arkr1998",
                                        "arkr1999",
                                        "arkr2004",
                                        "arkr2010",
                                        "arkr2013",
                                        "arkr2015"],

                    "Season Name"      : ["ARK 1995",
                                        "ARK 1996",
                                        "ARK 1997",
                                        "ARK 1998",
                                        "ARK 1999",
                                        "ARK 2004",
                                        "ARK 2010",
                                        "ARK 2013",
                                        "ARK 2015"],

                    "Radar System"     : "EMR",
                    "Radar Specs"    : "from ID",

                    "Platform"         : ["Polar 2",
                                        "Polar 2",
                                        "Polar 2",
                                        "Polar 2",
                                        "Polar 2",
                                        "Polar 2",
                                        "Polar 5",
                                        "Polar 5",
                                        "Polar 6"],

                    "Campaign Name"    : ["Nord GRIP",
                                        "Nord GRIP",
                                        "Nord GRIP",
                                        "NOGRAM",
                                        "NOGRAM",
                                        "North GRIP",
                                        "NEEM & EGIG",
                                        "Top 79 2013",
                                        "MABANG 1"],

                    "Campaign Name Long": ["Nord GRIP radar survey",
                                        "Nord GRIP radar survey",
                                        "Nord GRIP radar survey",
                                        "Northern Gravity, Radio Echo Sounding and Magnetics",
                                        "Northern Gravity, Radio Echo Sounding and Magnetics",
                                        "North GRIP",
                                        "North EEMian ice core project and International Glaciological Expedition to Greenland (EGIG) line",
                                        "opography and structure of the 79.5-Glacier and select areas in North-East Greenland",
                                        "Mass Balance of Glaciers in North Greenland (Part 1)"],

                    "Campaign PI"      : ["Daniel Steinhage (AWI)",
                                        "Daniel Steinhage (AWI)",
                                        "Daniel Steinhage (AWI)",
                                        "Heinrich Miller (AWI), Christoph Meyer (AWI)",
                                        "Heinrich Miller (AWI)",
                                        "Heinrich Miller (AWI)",
                                        "Heinrich Miller (AWI)",
                                        "Angelika Humbert (AWI)",
                                        "Angelika Humbert (AWI)"],

                    "Grant ID"          : ["None",
                                        "None",
                                        "None",
                                        "None",
                                        "None",
                                        "None",
                                        "AWI_PA_02017",
                                        "AWI_PA_02047",
                                        "AWI_PA_02072"],

                    "Pangaea Identifyer": ["Aeromag_1995",
                                        "NGRIP_1996",
                                        "NGRIP_1997",
                                        "NOGRAM_1998",
                                        "NOGRAM_1999",
                                        "NGRIP_2004",
                                        "NEEM_2010",
                                        "TOP79_5_2013",
                                        "MABANG_2015"],

                    "Radar Data DOI"    : ["None", #arkr1995
                                        "None", #arkr1996
                                        "None", #arkr1997
                                        "https://doi.org/10.1594/PANGAEA.974043", #arkr1998
                                        "https://doi.org/10.1594/PANGAEA.982470", #arkr1999
                                        "https://doi.org/10.1594/PANGAEA.982491", #arkr2004
                                        "https://doi.org/10.1594/PANGAEA.979248", #arkr2010
                                        "None", #arkr2013
                                        "https://doi.org/10.1594/PANGAEA.982577"  #arkr2015
                                        ]}

    ## ================================================================================================= ##
    ## ================================================================================================= ##


    list_UWB_ALL_dicts = []

    ###################
    # UWB ARK

    # ARK 2016 - Jakobshavn
    dict_UWB_arkr2016_Jakobshavn = {"Campaign Folder"    : "2016_Greenland_Polar6_rds_standard",
                                    "Season Name"        : "ARK 2016",
                                    "Radar System"       : "UWB",
                                    "Platform"           : "Polar 6",

                                    "Radar Specs"      : ["180-210 MHz bandwidth",
                                                            "180-210 MHz bandwidth",
                                                            "180-210 MHz bandwidth",
                                                            "180-210 MHz bandwidth",
                                                            "180-210 MHz bandwidth",
                                                            "180-210 MHz bandwidth",
                                                            "180-210 MHz bandwidth",
                                                            "180-210 MHz bandwidth",
                                                            "180-210 MHz bandwidth",
                                                            "180-210 MHz bandwidth",
                                                            "180-210 MHz bandwidth",
                                                            "180-210 MHz bandwidth",
                                                            "180-210 MHz bandwidth",
                                                            "150-520 MHz bandwidth",
                                                            "150-520 MHz bandwidth",
                                                            "180-210 MHz bandwidth",
                                                            "180-210 MHz bandwidth",
                                                            "180-210 MHz bandwidth",
                                                            "150-520 MHz bandwidth",
                                                            "180-210 MHz bandwidth",
                                                            "150-520 MHz bandwidth",
                                                            "150-520 MHz bandwidth",
                                                            "150-520 MHz bandwidth",
                                                            "150-520 MHz bandwidth",
                                                            "180-210 MHz bandwidth",
                                                            "180-210 MHz bandwidth",
                                                            "180-210 MHz bandwidth",
                                                            "150-520 MHz bandwidth"],

                                    "Segments"           : ["20160413_14",
                                                            "20160413_15",
                                                            "20160413_16",
                                                            "20160413_17",
                                                            "20160418_21",
                                                            "20160418_23",
                                                            "20160418_24",
                                                            "20160418_25",
                                                            "20160418_26",
                                                            "20160418_27",
                                                            "20160418_28",
                                                            "20160418_29",
                                                            "20160418_30",
                                                            "20160418_31",
                                                            "20160418_32",
                                                            "20160418_33",
                                                            "20160418_34",
                                                            "20160418_35",
                                                            "20160420_01",
                                                            "20160420_02",
                                                            "20160420_03",
                                                            "20160420_04",
                                                            "20160420_05",
                                                            "20160425_13",
                                                            "20160425_14",
                                                            "20160425_15",
                                                            "20160426_11",
                                                            "20160426_13"],

                                    "Pangaea Event"     : ["UWB_2016_1604130501",
                                                        "UWB_2016_1604130501",
                                                        "UWB_2016_1604130501",
                                                        "UWB_2016_1604130501",
                                                        "UWB_2016_1604180601",
                                                            "UWB_2016_1604180702",
                                                            "UWB_2016_1604180702",
                                                            "UWB_2016_1604180702",
                                                            "UWB_2016_1604180702",
                                                            "UWB_2016_1604180702",
                                                            "UWB_2016_1604180702",
                                                            "UWB_2016_1604180702",
                                                            "UWB_2016_1604180702",
                                                            "UWB_2016_1604180702",
                                                            "UWB_2016_1604180702",
                                                            "UWB_2016_1604180702",
                                                            "UWB_2016_1604180702",
                                                            "UWB_2016_1604180702",
                                                            "UWB_2016_1604200901",
                                                            "UWB_2016_1604200901",
                                                            "UWB_2016_1604200901",
                                                            "UWB_2016_1604200901",
                                                            "UWB_2016_1604200901",
                                                            "UWB_2016_1604251001",
                                                            "UWB_2016_1604251001",
                                                            "UWB_2016_1604251001",
                                                            "UWB_2016_1604261202",
                                                            "UWB_2016_1604261202"],

                                    "Campaign Name"      : "Jakobshavn Testflights",
                                    "Campaign Name Long" : "Testing of the new UWB and UWBM radars",
                                    "Campaign PI"        : "Heinrich Miller (AWI)",
                                    "Grant ID"           : "AWI_PA_02078",
                                    "Pangaea Identifyer" : "P6-201_UWB_2016",
                                    "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_arkr2016_Jakobshavn)

    # ARK 2016 - Hiawatha
    dict_UWB_arkr2016_HIAWATHA = {"Campaign Folder"    : "2016_Greenland_Polar6_rds_standard",
                                "Season Name"        : "ARK 2016",
                                "Radar System"       : "UWB",
                                "Platform"           : "Polar 6",
                                "Radar Specs"      : np.repeat("150-520 MHz bandwidth (24ch)", 11),

                                "Segments"           : ["20160512_02",
                                                        "20160512_03",
                                                        "20160512_04",
                                                        "20160512_05",
                                                        "20160512_07",
                                                        "20160516_02",
                                                        "20160516_04",
                                                        "20160517_03",
                                                        "20160517_04",
                                                        "20160517_06",
                                                        "20160517_08"],

                                "Pangaea Event"      : ["UWB_2016_1605121301",
                                                        "UWB_2016_1605121301",
                                                        "UWB_2016_1605121301",
                                                        "UWB_2016_1605121301",
                                                        "UWB_2016_1605121301",
                                                        "UWB_2016_1605161401",
                                                        "UWB_2016_1605161401",
                                                        "UWB_2016_1605171501",
                                                        "UWB_2016_1605171501",
                                                        "UWB_2016_1605171501",
                                                        "UWB_2016_1605171501"],

                                "Campaign Name"      : "Hiawatha",
                                "Campaign Name Long" : "The Hiawatha Structure – a New Young Impact Crater in Northern Greenland",
                                "Campaign PI"        : "Olaf Eisen (AWI), Joe MacGregor (NASA), Horst Machguth (Université de Fribourg)",
                                "Grant ID"           : "AWI_PA_02079",
                                "Pangaea Identifyer" : "P6-201_UWB_2016",
                                "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_arkr2016_HIAWATHA)

    # ARK 2016 - RESURV
    dict_UWB_arkr2016_RESURV79 = {"Campaign Folder"    : "2016_Greenland_Polar6_rds_standard",
                                "Season Name"        : "ARK 2016",
                                "Radar System"       : "UWB",
                                "Platform"           : "Polar 6",
                                "Radar Specs"      : np.repeat("150-520 MHz bandwidth", 4),

                                "Segments"           : ["20160620_01",
                                                        "20160627_08",
                                                        "20160629_02",
                                                        "20160629_03"],

                                "Pangaea Event"      : ["UWB_2016_1606201601",
                                                        "UWB_2016_1606271801",
                                                        "UWB_2016_1606291901",
                                                        "UWB_2016_1606291901"],

                                "Campaign Name"      : "RESURV 79",
                                "Campaign Name Long" : "Re-survey of the 79° Glacier (2016)",
                                "Campaign PI"        : "Angelika Humbert (AWI)",
                                "Grant ID"           : "AWI_PA_02080",
                                "Pangaea Identifyer" : "P6-201_UWB_2016",
                                "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_arkr2016_RESURV79)

    # ARK 2018 - RESURV79
    dict_UWB_arkr2018_RESURV79 = {"Campaign Folder"  : "2018_Greenland_Polar6_rds_standard",
                                "Season Name"        : "ARK 2018",
                                "Radar System"       : "UWB",
                                "Platform"           : "Polar 6",

                                "Radar Specs"        : ["150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth (swath mode)",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth (swath mode)",
                                                        "150-520 MHz bandwidth (swath mode)",
                                                        "180-210 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth (swath mode)",
                                                        "180-210 MHz bandwidth (swath mode)",
                                                        "180-210 MHz bandwidth (swath mode)",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth (swath mode)",
                                                        "180-210 MHz bandwidth (swath mode)",
                                                        "180-210 MHz bandwidth (swath mode)"],

                                "Segments"           : ["20180410_02",
                                                        "20180411_05",
                                                        "20180413_07",
                                                        "20180413_08",
                                                        "20180413_09",
                                                        "20180413_14",
                                                        "20180414_12",
                                                        "20180414_14",
                                                        "20180415_02",
                                                        "20180415_03",
                                                        "20180415_04",
                                                        "20180415_05",
                                                        "20180415_07",
                                                        "20180415_08",
                                                        "20180417_01",
                                                        "20180417_02",
                                                        "20180417_04",
                                                        "20180417_05",
                                                        "20180417_07",
                                                        "20180418_04",
                                                        "20180420_01",
                                                        "20180420_02",
                                                        "20180422_01",
                                                        "20180422_02",
                                                        "20180427_01",
                                                        "20180427_02",
                                                        "20180427_03"],

                                    "Pangaea Event"  : ["P6_211_RESURV79_2018_1804100301",
                                                        "P6_211_RESURV79_2018_1804110401",
                                                        "P6_211_RESURV79_2018_1804130501",
                                                        "P6_211_RESURV79_2018_1804130501",
                                                        "P6_211_RESURV79_2018_1804130501",
                                                        "P6_211_RESURV79_2018_1804130501",
                                                        "P6_211_RESURV79_2018_1804140601",
                                                        "P6_211_RESURV79_2018_1804140601",
                                                        "P6_211_RESURV79_2018_1804150701",
                                                        "P6_211_RESURV79_2018_1804150701",
                                                        "P6_211_RESURV79_2018_1804150701",
                                                        "P6_211_RESURV79_2018_1804150701",
                                                        "P6_211_RESURV79_2018_1804150701",
                                                        "P6_211_RESURV79_2018_1804150701",
                                                        "P6_211_RESURV79_2018_1804170801",
                                                        "P6_211_RESURV79_2018_1804170801",
                                                        "P6_211_RESURV79_2018_1804170801",
                                                        "P6_211_RESURV79_2018_1804170801",
                                                        "P6_211_RESURV79_2018_1804170801",
                                                        "P6_211_RESURV79_2018_1804180901",
                                                        "P6_211_RESURV79_2018_1804201001",
                                                        "P6_211_RESURV79_2018_1804201001",
                                                        "P6_211_RESURV79_2018_1804221101",
                                                        "P6_211_RESURV79_2018_1804221101",
                                                        "None",
                                                        "None",
                                                        "None"],

                                "Campaign Name"      : "RESURV 79",
                                "Campaign Name Long" : "Re-survey of the 79° Glacier (2018)",
                                "Campaign PI"        : "Angelika Humbert (AWI)",
                                "Grant ID"           : "AWI_PA_02095",
                                "Pangaea Identifyer" : "P6-211_RESURV79_2018",
                                "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_arkr2018_RESURV79)

    # ARK 2018 - FINEGIS
    dict_UWB_arkr2018_FINEGIS  = {"Campaign Folder"    : "2018_Greenland_Polar6_rds_standard",
                                "Season Name"        : "ARK 2018",
                                "Radar System"       : "UWB",
                                "Platform"           : "Polar 6",

                                "Radar Specs"      : ["150-520 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth"],

                                "Segments"           : ["20180414_08",
                                                        "20180414_09",
                                                        "20180415_06",
                                                        "20180418_01",
                                                        "20180418_03",
                                                        "20180423_03",
                                                        "20180423_04",
                                                        "20180423_05",
                                                        "20180423_06",
                                                        "20180423_07",
                                                        "20180423_08",
                                                        "20180423_09"],

                                "Pangaea Event"     : ["P6_211_RESURV79_2018_1804140601",
                                                        "P6_211_RESURV79_2018_1804140601",
                                                        "P6_211_RESURV79_2018_1804150701",
                                                        "P6_211_RESURV79_2018_1804180901",
                                                        "P6_211_RESURV79_2018_1804180901",
                                                        "P6_211_RESURV79_2018_1804231201",
                                                        "P6_211_RESURV79_2018_1804231201",
                                                        "P6_211_RESURV79_2018_1804231201",
                                                        "P6_211_RESURV79_2018_1804231201",
                                                        "P6_211_RESURV79_2018_1804231201",
                                                        "P6_211_RESURV79_2018_1804231201",
                                                        "P6_211_RESURV79_2018_1804231201"],

                                "Campaign Name"      : "FINEGIS",
                                "Campaign Name Long" : "Mapping Folds In the North-Eastern Greenland Ice Sheet",
                                "Campaign PI"        : "Daniela Jansen (AWI)",
                                "Grant ID"           : "AWI_PA_02096",
                                "Pangaea Identifyer" : "P6-211_RESURV79_2018",
                                "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_arkr2018_FINEGIS)

    # ARK 2018 - EGRIP-NOR-2018
    dict_UWB_arkr2018_EGRIP  = {"Campaign Folder"    : "2018_Greenland_Polar6_rds_standard",
                                "Season Name"        : "ARK 2018",
                                "Radar System"       : "UWB",
                                "Platform"           : "Polar 6",

                                "Radar Specs"      : ["180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth (swath mode)",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth (swath mode)",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth (swath mode)",
                                                        "180-210 MHz bandwidth"],

                                "Segments"           : ["20180508_01",
                                                        "20180508_02",
                                                        "20180508_03",
                                                        "20180508_04",
                                                        "20180508_05",
                                                        "20180508_06",
                                                        "20180509_01",
                                                        "20180510_01",
                                                        "20180510_02",
                                                        "20180511_01",
                                                        "20180512_01",
                                                        "20180512_02",
                                                        "20180514_01",
                                                        "20180514_03",
                                                        "20180515_01",
                                                        "20180517_01"],

                                "Pangaea Event"     : ["P6_211_EGRIP_NOR_2018_1805081501",
                                                        "P6_211_EGRIP_NOR_2018_1805081602",
                                                        "P6_211_EGRIP_NOR_2018_1805081602",
                                                        "P6_211_EGRIP_NOR_2018_1805081602",
                                                        "P6_211_EGRIP_NOR_2018_1805081602",
                                                        "P6_211_EGRIP_NOR_2018_1805081602",
                                                        "P6_211_EGRIP_NOR_2018_1805091701",
                                                        "P6_211_EGRIP_NOR_2018_1805101801",
                                                        "P6_211_EGRIP_NOR_2018_1805101902",
                                                        "P6_211_EGRIP_NOR_2018_1805112001",
                                                        "P6_211_EGRIP_NOR_2018_1805122101",
                                                        "P6_211_EGRIP_NOR_2018_1805122202",
                                                        "P6_211_EGRIP_NOR_2018_1805142301",
                                                        "P6_211_EGRIP_NOR_2018_1805142402",
                                                        "P6_211_EGRIP_NOR_2018_1805152501",
                                                        "P6_211_EGRIP_NOR_2018_1805172601"],

                                "Campaign Name"      : "EGRIP-NOR-2018",
                                "Campaign Name Long" : "Characterisation of the NEGIS Onset region with Radar around EastGRIP",
                                "Campaign PI"        : "Olaf Eisen (AWI)",
                                "Grant ID"           : "AWI_PA_02097",
                                "Pangaea Identifyer" : "P6-211_EGRIP_NOR_2018",
                                "Radar Data DOI"     : "https://doi.org/10.1594/PANGAEA.983191"}

    list_UWB_ALL_dicts.append(dict_UWB_arkr2018_EGRIP)

    # ARK 2021 - 79NG-EC
    dict_UWB_arkr2021_79NG_EC = {"Campaign Folder"   : "2021_Greenland_Polar5_rds_standard",
                                "Season Name"        : "ARK 2021",
                                "Radar System"       : "UWB",
                                "Platform"           : "Polar 5",

                                "Radar Specs"      : ["150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth"],

                                "Segments"           : ["20210727_01",
                                                        "20210727_02",
                                                        "20210727_03",
                                                        "20210727_07",
                                                        "20210727_09",
                                                        "20210728_01",
                                                        "20210728_03",
                                                        "20210729_01",
                                                        "20210729_02",
                                                        "20210730_01",
                                                        "20210730_04",
                                                        "20210730_05",
                                                        "20210803_01", 
                                                        "20210803_02",
                                                        "20210809_01", 
                                                        "20210809_03", 
                                                        "20210814_01"],

                                "Pangaea Event"     : ["P5_227_NG21_2021_2107270101",
                                                        "P5_227_NG21_2021_2107270101",
                                                        "P5_227_NG21_2021_2107270101",
                                                        "P5_227_NG21_2021_2107270101",
                                                        "P5_227_NG21_2021_2107270101",
                                                        "P5_227_NG21_2021_2107280201",
                                                        "P5_227_NG21_2021_2107280201",
                                                        "P5_227_NG21_2021_2107290301",
                                                        "P5_227_NG21_2021_2107290301",
                                                        "P5_227_NG21_2021_2107300401",
                                                        "P5_227_NG21_2021_2107300401",
                                                        "P5_227_NG21_2021_2107300401",
                                                        "P5_227_NG21_2021_2108030601", 
                                                        "P5_227_NG21_2021_2108030601", 
                                                        "P5_227_NG21_2021_2108090701", 
                                                        "P5_227_NG21_2021_2108090701"
                                                        "P5_227_NG21_2021_2108140801"],

                                "Campaign Name"      : "79NG-EC",
                                "Campaign Name Long" : "Englacial channels in 79°N Glacier, Greenland",
                                "Campaign PI"        : "Angelika Humbert (AWI)",
                                "Grant ID"           : "AWI_PA_02118",
                                "Pangaea Identifyer" : "P5-227_NG21_2021",
                                "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_arkr2021_79NG_EC)

    # ARK 2021 - HTRES
    dict_UWB_arkr2021_HTRES   = {"Campaign Folder"   : "2021_Greenland_Polar5_rds_standard",
                                "Season Name"        : "ARK 2021",
                                "Radar System"       : "UWB",
                                "Platform"           : "Polar 5",
                                "Radar Specs"      : ["150-520 MHz bandwidth"],
                                "Segments"           : ["20210731_01"],
                                "Pangaea Event"     : ["P5_227_NG21_2021_2107310501"],
                                "Campaign Name"      : "HTRES",
                                "Campaign Name Long" : "Hans Tausen Ice Cap Radio Echo Sounding",
                                "Campaign PI"        : "Angelika Humbert (AWI), Bo Vinther (University of Copenhagen)",
                                "Grant ID"           : "AWI_PA_02120",
                                "Pangaea Identifyer" : "P5-227_NG21_2021",
                                "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_arkr2021_HTRES)

    # ARK 2022 - NEGIS Folds
    dict_UWB_arkr2022_NEGSI_Folds_part1  = {"Campaign Folder"    : "2022_Greenland_Polar5_rds_standard",
                                            "Season Name"        : "ARK 2022",
                                            "Radar System"       : "UWB",
                                            "Platform"           : "Polar 5",
                                            "Radar Specs"      : ["180-210 MHz bandwidth"],
                                            "Segments"           : ["20220608_01"],
                                            "Pangaea Event"     : ["P5_233_NEGIS_2022_2206080401"],
                                            "Campaign Name"      : "NEGIS Folds",
                                            "Campaign Name Long" : "Folds adjacent to the North East Greenland Ice Stream",
                                            "Campaign PI"        : "Paul Bons (University of Tübingen)",
                                            "Grant ID"           : "AWI_PA_02130",
                                            "Pangaea Identifyer" : "P5-233_NEGIS_2022",
                                            "Radar Data DOI"     : "None"}
    list_UWB_ALL_dicts.append(dict_UWB_arkr2022_NEGSI_Folds_part1)

    dict_UWB_arkr2022_NEGSI_Folds_part2  = {"Campaign Folder"    : "2022_Greenland_Polar5_rds_standard_VV",
                                            "Season Name"        : "ARK 2022",
                                            "Radar System"       : "UWB",
                                            "Platform"           : "Polar 5",
                                            "Radar Specs"      : ["180-210 MHz bandwidth (4ch)"],
                                            "Segments"           : ["20220616_01"],
                                            "Pangaea Event"     : ["P5_233_NEGIS_2022_2206161001"],
                                            "Campaign Name"      : "NEGIS Folds",
                                            "Campaign Name Long" : "Folds adjacent to the North East Greenland Ice Stream",
                                            "Campaign PI"        : "Paul Bons (University of Tübingen)",
                                            "Grant ID"           : "AWI_PA_02130",
                                            "Pangaea Identifyer" : "P5-233_NEGIS_2022",
                                            "Radar Data DOI"     : "None"}
    list_UWB_ALL_dicts.append(dict_UWB_arkr2022_NEGSI_Folds_part2)

    # ARK 2022 - NEGIS-FLOW
    dict_UWB_arkr2022_NEGSI_FLOW     = {"Campaign Folder"    : "2022_Greenland_Polar5_rds_standard",
                                        "Season Name"        : "ARK 2022",
                                        "Radar System"       : "UWB",
                                        "Platform"           : "Polar 5",
                                        
                                        "Radar Specs"      : ["180-210 MHz bandwidth",
                                                                "180-210 MHz bandwidth",
                                                                "180-210 MHz bandwidth",
                                                                "180-210 MHz bandwidth",
                                                                "180-210 MHz bandwidth (swath mode)",
                                                                "180-210 MHz bandwidth",
                                                                "180-210 MHz bandwidth",
                                                                "180-210 MHz bandwidth"],

                                        "Segments"           : ["20220603_02",
                                                                "20220607_01",
                                                                "20220609_01",
                                                                "20220610_01",
                                                                "20220610_02",
                                                                "20220611_01",
                                                                "20220612_01",
                                                                "20220612_02"],

                                        "Pangaea Event"     : ["P5_233_NEGIS_2022_2206030202",
                                                                "P5_233_NEGIS_2022_2206070301",
                                                                "P5_233_NEGIS_2022_2206090501",
                                                                "P5_233_NEGIS_2022_2206100601",
                                                                "P5_233_NEGIS_2022_2206100601",
                                                                "P5_233_NEGIS_2022_2206110701",
                                                                "P5_233_NEGIS_2022_2206120801",
                                                                "P5_233_NEGIS_2022_2206120801"],
                                        
                                        "Campaign Name"      : "NEGIS-FLOW",
                                        "Campaign Name Long" : "Flow characteristics of the North East Greenland Ice Stream",
                                        "Campaign PI"        : "Olaf Eisen (AWI)",
                                        "Grant ID"           : "AWI_PA_02128",
                                        "Pangaea Identifyer" : "P5-233_NEGIS_2022",
                                        "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_arkr2022_NEGSI_FLOW)

    # ARK 2022 - NEGIS-ANISO
    dict_UWB_arkr2022_NEGSI_ANISO     = {"Campaign Folder"   : "2022_Greenland_Polar5_rds_standard_VV",
                                        "Season Name"        : "ARK 2022",
                                        "Radar System"       : "UWB",
                                        "Platform"           : "Polar 5",
                                        
                                        "Radar Specs"        : ["180-210 MHz bandwidth (polarized mode)"
                                                                "180-210 MHz bandwidth (polarized mode)",
                                                                "180-210 MHz bandwidth (polarized mode)",
                                                                "180-210 MHz bandwidth (polarized mode)",
                                                                "180-210 MHz bandwidth (polarized mode)",
                                                                "180-210 MHz bandwidth (polarized mode)",
                                                                "180-210 MHz bandwidth (polarized mode)"],

                                        "Segments"           : ["20220609_02",
                                                                "20220612_03",
                                                                "20220613_01",
                                                                "20220613_02",
                                                                "20220613_03",
                                                                "20220613_04",
                                                                "20220613_05"],

                                        "Pangaea Event"      : ["P5_233_NEGIS_2022_2206090501",
                                                                "P5_233_NEGIS_2022_2206120801",
                                                                "P5_233_NEGIS_2022_2206130901",
                                                                "P5_233_NEGIS_2022_2206130901",
                                                                "P5_233_NEGIS_2022_2206130901",
                                                                "P5_233_NEGIS_2022_2206130901",
                                                                "P5_233_NEGIS_2022_2206130901"],
                                        
                                        "Campaign Name"      : "NEGIS-ANISO",
                                        "Campaign Name Long" : "Anisotropy of the North-East Greenland ice stream",
                                        "Campaign PI"        : "Olaf Eisen (AWI)",
                                        "Grant ID"           : "AWI_PA_02129",
                                        "Pangaea Identifyer" : "P5-233_NEGIS_2022",
                                        "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_arkr2022_NEGSI_ANISO)

    # ARK 2023 - Müller Ice Cap
    dict_UWB_arkr2023_MIC = {"Campaign Folder"    : "2023_Canada_Polar5_rds_standard",
                            "Season Name"        : "ARK 2023",
                            "Radar System"       : "UWB",
                            "Platform"           : "Polar 5",

                            "Radar Specs"      : ["150-520 MHz bandwidth",
                                                    "150-520 MHz bandwidth",
                                                    "180-210 MHz bandwidth"],

                            "Segments"           : ["20230508_09",
                                                    "20230508_10",
                                                    "20230509_03"],

                            "Pangaea Event"     : ["P5_240_MullerIceCap_2023_202305080101",
                                                    "P5_240_MullerIceCap_2023_202305080101",
                                                    "P5_240_MullerIceCap_2023_202305090201"],

                            "Campaign Name"      : "Müller Ice Cap",
                            "Campaign Name Long" : "Selection of a drill site on the Müller Ice Cap, Axel Heibergs Island, Canada",
                            "Campaign PI"        : "Daniel Steinhage (AWI), Dorthe Dahl-Jensen (University of Manitoba)",
                            "Grant ID"           : "AWI_PA_02140",
                            "Pangaea Identifyer" : "P5-240_MullerIceCap_2023",
                            "Radar Data DOI"     : "https://doi.org/10.1594/PANGAEA.980473"}

    list_UWB_ALL_dicts.append(dict_UWB_arkr2023_MIC)

    # ARK 2024 - SLOGIS
    dict_UWB_arkr2024_SLOGIS     = {"Campaign Folder"    : "2024_Greenland_Polar6_rds_standard",
                                    "Season Name"        : "ARK 2024",
                                    "Radar System"       : "UWB",
                                    "Platform"           : "Polar 6",
                                    "Radar Specs"      : np.repeat("150-520 MHz bandwidth", 11),

                                    "Segments"           : ["20240727_01",
                                                            "20240727_03",
                                                            "20240731_01",
                                                            "20240801_01",
                                                            "20240802_01",
                                                            "20240802_02",
                                                            "20240803_01",
                                                            "20240803_02",
                                                            "20240803_03",
                                                            "20240803_04",
                                                            "20240803_05"],

                                    "Pangaea Event"     : ["P6-252_SLOGIS_2024_2407270101",
                                                            "P6-252_SLOGIS_2024_2407270101",
                                                            "P6-252_SLOGIS_2024_2407310201",
                                                            "P6-252_SLOGIS_2024_2408010301",
                                                            "P6-252_SLOGIS_2024_2408020401",
                                                            "P6-252_SLOGIS_2024_2408020401",
                                                            "P6-252_SLOGIS_2024_2408030501",
                                                            "P6-252_SLOGIS_2024_2408030501",
                                                            "P6-252_SLOGIS_2024_2408030501",
                                                            "P6-252_SLOGIS_2024_2408030501",
                                                            "P6-252_SLOGIS_2024_2408030501"],

                                    "Campaign Name"      : "SLOGIS",
                                    "Campaign Name Long" : "Supraglacial Lakes On the Greenland Ice Sheet – 79°N Glacier",
                                    "Campaign PI"        : "Angelika Humbert (AWI)",
                                    "Grant ID"           : "AWI_PA_02154",
                                    "Pangaea Identifyer" : "P6-252_SLOGIS_2024",
                                    "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_arkr2024_SLOGIS)

    # ARK 2024 - Polarmonitor
    dict_UWB_arkr2024_Polarmonitor   = {"Campaign Folder"    : "2024_Greenland_Polar6_rds_standard",
                                        "Season Name"        : "ARK 2024",
                                        "Radar System"       : "UWB",
                                        "Platform"           : "Polar 6",
                                        "Radar Specs"      : np.repeat("150-520 MHz bandwidth", 12),

                                        "Segments"           : ["20240817_01",
                                                                "20240819_01",
                                                                "20240819_02",
                                                                "20240819_03",
                                                                "20240820_01",
                                                                "20240821_01",
                                                                "20240821_02",
                                                                "20240822_01",
                                                                "20240828_01"],

                                        "Pangaea Event"     : ["P6-253_Polarmonitor_2024_2408170101",
                                                                "P6-253_Polarmonitor_2024_2408190201",
                                                                "P6-253_Polarmonitor_2024_2408190201",
                                                                "P6-253_Polarmonitor_2024_2408190201",
                                                                "P6-253_Polarmonitor_2024_2408200301",
                                                                "P6-253_Polarmonitor_2024_2408210401",
                                                                "P6-253_Polarmonitor_2024_2408210502",
                                                                "P6-253_Polarmonitor_2024_2408220601",
                                                                "P6-253_Polarmonitor_2024_2408280701"],

                                        "Campaign Name"      : "Polarmonitor",
                                        "Campaign Name Long" : "Polarmonitor (2024)",
                                        "Campaign PI"        : "Angelika Humbert (AWI)",
                                        "Grant ID"           : "AWI_PA_02155",
                                        "Pangaea Identifyer" : "P6-253_Polarmonitor_2024",
                                        "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_arkr2024_Polarmonitor)

    # ARK 2024 - ATIWIK
    dict_UWB_arkr2024_ATIWIK         = {"Campaign Folder"    : "2024_Greenland_Polar6_rds_standard",
                                        "Season Name"        : "ARK 2024",
                                        "Radar System"       : "UWB",
                                        "Platform"           : "Polar 6",
                                        "Radar Specs"      : np.repeat("150-520 MHz bandwidth", 3),

                                        "Segments"           : ["20240901_01",
                                                                "20240902_01",
                                                                "20240902_02"],

                                        "Pangaea Event"     : ["P6-254_ATIWIK_2024_2409010102",
                                                                "P6-254_ATIWIK_2024_2409020201",
                                                                "P6-254_ATIWIK_2024_2409020201"],

                                        "Campaign Name"      : "ATIWIK",
                                        "Campaign Name Long" : "Arctic tidal wetlands in Kalaallit Nunaat (Greenland): Geomorphological properties and ecological values of polar sedimentary shorelines",
                                        "Campaign PI"        : "Angelika Humbert (AWI), Lasse Sander (AWI)",
                                        "Grant ID"           : "AWI_PA_02156",
                                        "Pangaea Identifyer" : "P6-254_ATIWIK_2024",
                                        "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_arkr2024_ATIWIK)


    # ARK 2024 - PERMA-X Alaska 2024
    dict_UWB_arkr2024_PERMA_X         = {"Campaign Folder"   : "2024_Alaska_Polar6_rds_standard",
                                        "Season Name"        : "ARK 2024",
                                        "Radar System"       : "UWB",
                                        "Platform"           : "Polar 6",
                                        "Radar Specs"        : np.repeat("150-520 MHz bandwidth", 15),

                                        "Segments"           : ["20240716_02",
                                                                "20240716_03",
                                                                "20240716_04",
                                                                "20240716_05",
                                                                "20240716_06",
                                                                "20240716_07",
                                                                "20240716_08",
                                                                "20240716_09",
                                                                "20240716_10",
                                                                "20240716_11",
                                                                "20240718_01",
                                                                "20240718_02",
                                                                "20240718_03",
                                                                "20240718_04",
                                                                "20240718_05"],

                                        "Pangaea Event"      : ["P6-249_Perma-X_2024_2407161401",
                                                                "P6-249_Perma-X_2024_2407161401",
                                                                "P6-249_Perma-X_2024_2407161401",
                                                                "P6-249_Perma-X_2024_2407161401",
                                                                "P6-249_Perma-X_2024_2407161401",
                                                                "P6-249_Perma-X_2024_2407161401",
                                                                "P6-249_Perma-X_2024_2407161401",
                                                                "P6-249_Perma-X_2024_2407161401",
                                                                "P6-249_Perma-X_2024_2407161401",
                                                                "P6-249_Perma-X_2024_2407161401",
                                                                "P6-249_Perma-X_2024_2407171501",
                                                                "P6-249_Perma-X_2024_2407171501",
                                                                "P6-249_Perma-X_2024_2407171501",
                                                                "P6-249_Perma-X_2024_2407171501",
                                                                "P6-249_Perma-X_2024_2407171501"],

                                        "Campaign Name"      : "PERMA-X Alaska 2024",
                                        "Campaign Name Long" : "Permafrost thaw dynamics, coastal erosion, vegetation and carbon fluxes on the North Slope of Alaska and in Western Boreal to Subarctic Canada",
                                        "Campaign PI"        : "Guido Grosse (AWI)",
                                        "Grant ID"           : "AWI_PA_02151",
                                        "Pangaea Identifyer" : "P6-249_Perma-X_2024",
                                        "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_arkr2024_PERMA_X)

    ###############
    # UWB ANT

    # ANT 2017 - MABAJU 2
    dict_UWB_antr2017_MABAJU2       =  {"Campaign Folder"    : "2017_Antarctica_Polar6_rds_standard",
                                        "Season Name"        : "ANT 2016/17",
                                        "Radar System"       : "UWB",
                                        "Platform"           : "Polar 6",
                                        "Radar Specs"      : np.repeat("150-520 MHz bandwidth", 33),

                                        "Segments"           : ["20170118_02",
                                                                "20170118_03",
                                                                "20170118_04",
                                                                "20170118_05",
                                                                "20170118_10",
                                                                "20170118_11",
                                                                "20170118_12",
                                                                "20170118_13",
                                                                "20170118_14",
                                                                "20170121_08",
                                                                "20170122_02",
                                                                "20170122_03",
                                                                "20170122_04",
                                                                "20170122_05",
                                                                "20170122_06",
                                                                "20170122_07",
                                                                "20170122_08",
                                                                "20170122_09",
                                                                "20170122_10",
                                                                "20170122_11",
                                                                "20170122_12",
                                                                "20170122_13",
                                                                "20170123_02",
                                                                "20170123_03",
                                                                "20170123_04",
                                                                "20170123_05",
                                                                "20170123_06",
                                                                "20170124_02",
                                                                "20170124_03",
                                                                "20170124_04",
                                                                "20170124_05",
                                                                "20170124_06",
                                                                "20170124_07"],

                                        "Pangaea Event"     : ["P6_204_ANT_UWB_2016_2017_1701180201",
                                                               "P6_204_ANT_UWB_2016_2017_1701180201",
                                                               "P6_204_ANT_UWB_2016_2017_1701180201",
                                                               "P6_204_ANT_UWB_2016_2017_1701180201",
                                                               "P6_204_ANT_UWB_2016_2017_1701180201",
                                                               "P6_204_ANT_UWB_2016_2017_1701180201",
                                                               "P6_204_ANT_UWB_2016_2017_1701180201",
                                                               "P6_204_ANT_UWB_2016_2017_1701180201",
                                                               "P6_204_ANT_UWB_2016_2017_1701180201",
                                                               "P6_204_ANT_UWB_2016_2017_1701210301",
                                                               "P6_204_ANT_UWB_2016_2017_1701220401",
                                                                "P6_204_ANT_UWB_2016_2017_1701220401",
                                                                "P6_204_ANT_UWB_2016_2017_1701220401",
                                                                "P6_204_ANT_UWB_2016_2017_1701220401",
                                                                "P6_204_ANT_UWB_2016_2017_1701220502",
                                                                "P6_204_ANT_UWB_2016_2017_1701220502",
                                                                "P6_204_ANT_UWB_2016_2017_1701220502",
                                                                "P6_204_ANT_UWB_2016_2017_1701220502",
                                                                "P6_204_ANT_UWB_2016_2017_1701220603",
                                                                "P6_204_ANT_UWB_2016_2017_1701220603",
                                                                "P6_204_ANT_UWB_2016_2017_1701220603",
                                                                "P6_204_ANT_UWB_2016_2017_1701220603",
                                                                "P6_204_ANT_UWB_2016_2017_1701230701",
                                                                "P6_204_ANT_UWB_2016_2017_1701230701",
                                                                "P6_204_ANT_UWB_2016_2017_1701230701",
                                                                "P6_204_ANT_UWB_2016_2017_1701230701",
                                                                "P6_204_ANT_UWB_2016_2017_1701230701",
                                                                "P6_204_ANT_UWB_2016_2017_1701240801",
                                                                "P6_204_ANT_UWB_2016_2017_1701240801",
                                                                "P6_204_ANT_UWB_2016_2017_1701240801",
                                                                "P6_204_ANT_UWB_2016_2017_1701240801",
                                                                "P6_204_ANT_UWB_2016_2017_1701240801",
                                                                "P6_204_ANT_UWB_2016_2017_1701240801"],

                                        "Campaign Name"      : "MABAJU 2",
                                        "Campaign Name Long" : "Mass balance of Jutulstraumen (Part 2)",
                                        "Campaign PI"        : "Angelika Humbert (AWI)",
                                        "Grant ID"           : "AWI_PA_02085",
                                        "Pangaea Identifyer" : "P6-204_ANT_UWB_2016_2017",
                                        "Radar Data DOI"     : "https://doi.org/10.1594/PANGAEA.982758"}

    list_UWB_ALL_dicts.append(dict_UWB_antr2017_MABAJU2)

    # ANT 2018/19 - JURAS
    dict_UWB_antr2019_JURAS  = {"Campaign Folder"    : "2019_Antarctica_Polar6_rds_standard",
                                "Season Name"        : "ANT 2018/19",
                                "Radar System"       : "UWB",
                                "Radar Specs"      : np.repeat("180-210 MHz bandwidth", 19),
                                "Platform"           : "Polar 6",

                                "Segments"           : ["20181220_05",
                                                        "20181220_06",
                                                        "20181220_07",
                                                        "20181222_01",
                                                        "20181226_01",
                                                        "20181226_02",
                                                        "20181226_03",
                                                        "20181227_01",
                                                        "20181227_02",
                                                        "20181227_03",
                                                        "20181227_04",
                                                        "20181230_01",
                                                        "20181230_02",
                                                        "20181230_08",
                                                        "20181230_09",
                                                        "20190105_01",
                                                        "20190105_02",
                                                        "20190105_03",
                                                        "20190105_04"],

                                "Pangaea Event"      : ["P6_215_UWB_2018_1812220201",
                                                        "P6_215_UWB_2018_1812220201",
                                                        "P6_215_UWB_2018_1812220201",
                                                        "P6_215_UWB_2018_1812220301",
                                                        "P6_215_UWB_2018_1812260501",
                                                        "P6_215_UWB_2018_1812260501",
                                                        "P6_215_UWB_2018_1812260602",
                                                        "P6_215_UWB_2018_1812270701",
                                                        "P6_215_UWB_2018_1812270701",
                                                        "P6_215_UWB_2018_1812270802",
                                                        "P6_215_UWB_2018_1812270802",
                                                        "P6_215_UWB_2018_1812300901",
                                                        "P6_215_UWB_2018_1812300901",
                                                        "P6_215_UWB_2018_1812300901",
                                                        "P6_215_UWB_2018_1812300901",
                                                        "P6_215_UWB_2018_1901051101",
                                                        "P6_215_UWB_2018_1901051101",
                                                        "P6_215_UWB_2018_1901051101",
                                                        "P6_215_UWB_2018_1901051101"],

                                "Campaign Name"      : "JuRaS",
                                "Campaign Name Long" : "Jutulstraumen Radar Stratigraphy",
                                "Campaign PI"        : "Daniela Jansen (AWI)",
                                "Grant ID"           : "AWI_PA_02105",
                                "Pangaea Identifyer" : "P6-215_UWB_2018",
                                "Radar Data DOI"     : "https://doi.org/10.1594/PANGAEA.976759"}

    list_UWB_ALL_dicts.append(dict_UWB_antr2019_JURAS)

    # ANT 2018/19 - CHIRP
    dict_UWB_antr2019_CHIRP  = {"Campaign Folder"    : "2019_Antarctica_Polar6_rds_standard",
                                "Season Name"        : "ANT 2018/19",
                                "Radar System"       : "UWB",
                                "Radar Specs"      : np.repeat("150-520 MHz bandwidth", 5),
                                "Platform"           : "Polar 6",

                                "Segments"           : ["20181231_01",
                                                        "20190106_01",
                                                        "20190106_02",
                                                        "20190107_01",
                                                        "20190107_02"],

                                "Pangaea Event"     : ["P6_215_UWB_2018_1812311001",
                                                        "P6_215_UWB_2018_1901061201",
                                                        "P6_215_UWB_2018_1901061302",
                                                        "P6_215_UWB_2018_1901071401",
                                                        "P6_215_UWB_2018_1901071401"],

                                "Campaign Name"      : "CHIRP",
                                "Campaign Name Long" : "Channel and Ice Rise Project in DML",
                                "Campaign PI"        : "Olaf Eisen (AWI), Reinhard Drews (University of Tübingen)",
                                "Grant ID"           : "AWI_PA_02106",
                                "Pangaea Identifyer" : "P6-215_UWB_2018",
                                "Radar Data DOI"     : "https://doi.org/10.1594/PANGAEA.982515"}

    list_UWB_ALL_dicts.append(dict_UWB_antr2019_CHIRP)

    # ANT 2023/24 - CHARISO
    dict_UWB_antr2024_CHARISO = {"Campaign Folder"    : "2024_Antarctica_Polar6_rds_standard",
                                "Season Name"        : "ANT 2023/24",
                                "Radar System"       : "UWB",
                                "Radar Specs"      : np.repeat("180-210 MHz bandwidth", 3),
                                "Platform"           : "Polar 6",

                                "Segments"           : ["20231203_04",
                                                        "20231204_01",
                                                        "20231204_03"],

                                "Pangaea Event"     : ["P6_244_ANT_23_24_2312031001",
                                                        "P6_244_ANT_23_24_2312041201",
                                                        "P6_244_ANT_23_24_2312041201"],

                                "Campaign Name"      : "CHARISO",
                                "Campaign Name Long" : "Isochrones in East Antarctic outflow regions (Swiss Antarctic Research Initiative)",
                                "Campaign PI"        : "Olaf Eisen (AWI), Johannes Sutter (University of Bern)",
                                "Grant ID"           : "AWI_PA_02145",
                                "Pangaea Identifyer" : "P6-244_ANT_23_24",
                                "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_antr2024_CHARISO)

    # ANT 2023/24 - SNaCC DML
    dict_UWB_antr2024_SNACC_DML = {"Campaign Folder"    : "2024_Antarctica_Polar6_rds_standard",
                                "Season Name"        : "ANT 2023/24",
                                "Radar System"       : "UWB",
                                "Radar Specs"      : np.repeat("150-520 MHz bandwidth", 3),
                                "Platform"           : "Polar 6",

                                "Segments"           : ["20231129_01",
                                                        "20231211_01",
                                                        "20231211_02"],

                                "Pangaea Event"     : ["P6_244_ANT_23_24_2311290701",
                                                        "P6_244_ANT_23_24_2312111501",
                                                        "P6_244_ANT_23_24_2312111602"],

                                "Campaign Name"      : "SNaCC DML",
                                "Campaign Name Long" : "Dronning Maud Land Snow Accumulation rate reconstructions",
                                "Campaign PI"        : "Alexandra Zuhr (University of Tübingen)",
                                "Grant ID"           : "AWI_PA_02146",
                                "Pangaea Identifyer" : "P6-244_ANT_23_24",
                                "Radar Data DOI"     : "https://doi.pangaea.de/10.1594/PANGAEA.981128"}

    list_UWB_ALL_dicts.append(dict_UWB_antr2024_SNACC_DML)

    # ANT 2023/24 - RINGS DML
    dict_UWB_antr2024_RINGS_DML = {"Campaign Folder"    : "2024_Antarctica_Polar6_rds_standard",
                                "Season Name"        : "ANT 2023/24",
                                "Radar System"       : "UWB",
                                "Platform"           : "Polar 6",

                                "Radar Specs"      : ["150-520 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth (swath mode)"],
                                
                                "Segments"           : ["20231125_15",
                                                        "20231126_01",
                                                        "20231126_02",
                                                        "20231126_04",
                                                        "20231126_05",
                                                        "20231127_01",
                                                        "20231128_03",
                                                        "20231128_04",
                                                        "20231130_01",
                                                        "20231130_02",
                                                        "20231130_03",
                                                        "20231130_05",
                                                        "20231203_05",
                                                        "20231204_04",
                                                        "20231204_05"],

                                "Pangaea Event"     : ["P6_244_ANT_23_24_2311250202",
                                                        "P6_244_ANT_23_24_2311260301",
                                                        "P6_244_ANT_23_24_2311260301",
                                                        "P6_244_ANT_23_24_2311260301",
                                                        "P6_244_ANT_23_24_2311260301",
                                                        "P6_244_ANT_23_24_2311270401",
                                                        "P6_244_ANT_23_24_2311280602",
                                                        "P6_244_ANT_23_24_2311280602",
                                                        "P6_244_ANT_23_24_2311300801",
                                                        "P6_244_ANT_23_24_2311300801",
                                                        "P6_244_ANT_23_24_2311300801",
                                                        "P6_244_ANT_23_24_2311300801",
                                                        "P6_244_ANT_23_24_2312031102",
                                                        "P6_244_ANT_23_24_2312041302",
                                                        "P6_244_ANT_23_24_2312041302"],

                                "Campaign Name"      : "RINGS DML",
                                "Campaign Name Long" : "Mapping ice thickness near the grounding line in Dronning Maud Land for SCAR-RINGS",
                                "Campaign PI"        : "Olaf Eisen (AWI), Kenichi Matsuoka (NPI)",
                                "Grant ID"           : "AWI_PA_02147",
                                "Pangaea Identifyer" : "P6-244_ANT_23_24",
                                "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_antr2024_RINGS_DML)

    # ANT 2023/24 - CHIRP 2
    dict_UWB_antr2024_CHIRP2 = {"Campaign Folder"    : "2024_Antarctica_Polar6_rds_standard",
                                "Season Name"        : "ANT 2023/24",
                                "Radar System"       : "UWB",
                                "Platform"           : "Polar 6",

                                "Radar Specs"      : ["150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "180-210 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "150-520 MHz bandwidth",
                                                        "180-210 MHz bandwidth"],

                                "Segments"           : ["20231125_12",
                                                        "20231125_13",
                                                        "20231125_14",
                                                        "20231128_01",
                                                        "20231128_02",
                                                        "20231130_04",
                                                        "20231202_01",
                                                        "20231202_02",
                                                        "20231204_02",
                                                        "20231210_01"],

                                "Pangaea Event"     : ["P6_244_ANT_23_24_2311250101",
                                                        "P6_244_ANT_23_24_2311250101",
                                                        "P6_244_ANT_23_24_2311250101",
                                                        "P6_244_ANT_23_24_2311280501",
                                                        "P6_244_ANT_23_24_2311280501",
                                                        "P6_244_ANT_23_24_2311300801",
                                                        "P6_244_ANT_23_24_2312020901",
                                                        "P6_244_ANT_23_24_2312020901",
                                                        "P6_244_ANT_23_24_2312041201",
                                                        "P6_244_ANT_23_24_2312101401"],

                                "Campaign Name"      : "CHIRP 2",
                                "Campaign Name Long" : "Channel and Ice Rise Project in DML 2",
                                "Campaign PI"        : "Olaf Eisen (AWI), Reinhard Drews (University of Tübingen)",
                                "Grant ID"           : "None",
                                "Pangaea Identifyer" : "P6-244_ANT_23_24",
                                "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_antr2024_CHIRP2)


    # ANT 2024/25 - RINGS
    dict_UWB_antr2025_RINGS = {"Campaign Folder"     : "2025_Antarctica_Polar6_rds_standard",
                                "Season Name"        : "ANT 2024/25",
                                "Radar System"       : "UWB",
                                "Platform"           : "Polar 6",
                                "Radar Specs"      : np.repeat("150-520 MHz bandwidth", 54),

                                "Segments"           : ["20241204_01",
                                                        "20241204_03",
                                                        "20241204_04",
                                                        "20241204_05",
                                                        "20241205_01",
                                                        "20241205_02",
                                                        "20241205_03",
                                                        "20241205_04",
                                                        "20241205_05",
                                                        "20241206_01",
                                                        "20241206_02",
                                                        "20241206_03",
                                                        "20241206_04",
                                                        "20241206_05",
                                                        "20241206_06",
                                                        "20241207_01",
                                                        "20241207_02",
                                                        "20241207_03",
                                                        "20241207_04",
                                                        "20241208_01",
                                                        "20241208_02",
                                                        "20241208_03",
                                                        "20241208_04",
                                                        "20241210_01",
                                                        "20241210_02",
                                                        "20241210_03",
                                                        "20241212_01",
                                                        "20241212_04",
                                                        "20241212_05",
                                                        "20241212_06",
                                                        "20241212_08",
                                                        "20241213_01",
                                                        "20241215_01",
                                                        "20241215_02",
                                                        "20241215_03",
                                                        "20241216_01",
                                                        "20241217_01",
                                                        "20241218_01",
                                                        "20241218_02",
                                                        "20241218_03",
                                                        "20241221_01",
                                                        "20241221_02",
                                                        "20241221_03",
                                                        "20241222_01",
                                                        "20241222_02",
                                                        "20250101_01",
                                                        "20250101_02",
                                                        "20250101_03",
                                                        "20250105_02",
                                                        "20250106_01",
                                                        "20250107_01",
                                                        "20250107_03",
                                                        "20250107_04",
                                                        "20250107_05"],

                                "Pangaea Event"     : ["P6-255_ANT_2024_2025_2412040101",
                                                        "P6-255_ANT_2024_2025_2412040101",
                                                        "P6-255_ANT_2024_2025_2412040101",
                                                        "P6-255_ANT_2024_2025_2412040101",
                                                        "P6-255_ANT_2024_2025_2412050201",
                                                        "P6-255_ANT_2024_2025_2412050201",
                                                        "P6-255_ANT_2024_2025_2412050201",
                                                        "P6-255_ANT_2024_2025_2412050201",
                                                        "P6-255_ANT_2024_2025_2412050201",
                                                        "P6-255_ANT_2024_2025_2412060301",
                                                        "P6-255_ANT_2024_2025_2412060301",
                                                        "P6-255_ANT_2024_2025_2412060301",
                                                        "P6-255_ANT_2024_2025_2412060301",
                                                        "P6-255_ANT_2024_2025_2412060301",
                                                        "P6-255_ANT_2024_2025_2412060301",
                                                        "P6-255_ANT_2024_2025_2412070401",
                                                        "P6-255_ANT_2024_2025_2412070401",
                                                        "P6-255_ANT_2024_2025_2412070401",
                                                        "P6-255_ANT_2024_2025_2412070401",
                                                        "P6-255_ANT_2024_2025_2412080501",
                                                        "P6-255_ANT_2024_2025_2412080501",
                                                        "P6-255_ANT_2024_2025_2412080501",
                                                        "P6-255_ANT_2024_2025_2412080501",
                                                        "P6-255_ANT_2024_2025_2412100601",
                                                        "P6-255_ANT_2024_2025_2412100702",
                                                        "P6-255_ANT_2024_2025_2412100702",
                                                        "P6-255_ANT_2024_2025_2412120801",
                                                        "P6-255_ANT_2024_2025_2412120801",
                                                        "P6-255_ANT_2024_2025_2412120801",
                                                        "P6-255_ANT_2024_2025_2412120801",
                                                        "P6-255_ANT_2024_2025_2412120801",
                                                        "P6-255_ANT_2024_2025_2412130901",
                                                        "P6-255_ANT_2024_2025_2412151101",
                                                        "P6-255_ANT_2024_2025_2412151101",
                                                        "P6-255_ANT_2024_2025_2412151101",
                                                        "P6-255_ANT_2024_2025_2412161201",
                                                        "P6-255_ANT_2024_2025_2412171301",
                                                        "P6-255_ANT_2024_2025_2412181401",
                                                        "P6-255_ANT_2024_2025_2412181401",
                                                        "P6-255_ANT_2024_2025_2412181401",
                                                        "P6-255_ANT_2024_2025_2412211701",
                                                        "P6-255_ANT_2024_2025_2412211701",
                                                        "P6-255_ANT_2024_2025_2412211701",
                                                        "P6-255_ANT_2024_2025_2412221801",
                                                        "P6-255_ANT_2024_2025_2412221801",
                                                        "P6-255_ANT_2024_2025_2501012002",
                                                        "P6-255_ANT_2024_2025_2501012002",
                                                        "P6-255_ANT_2024_2025_2501012002",
                                                        "P6-255_ANT_2024_2025_2501052201",
                                                        "P6-255_ANT_2024_2025_2501062301",
                                                        "P6-255_ANT_2024_2025_2501072401",
                                                        "P6-255_ANT_2024_2025_2501072401",
                                                        "P6-255_ANT_2024_2025_2501072502",
                                                        "P6-255_ANT_2024_2025_2501072502"],

                                "Campaign Name"      : "RINGS EL",
                                "Campaign Name Long" : "Mapping ice thickness near the grounding line in Dronning Maud Land for SCAR-RINGS",
                                "Campaign PI"        : "Olaf Eisen (AWI), Kenichi Matsuoka (NPI)",
                                "Grant ID"           : "AWI_PA_02157",
                                "Pangaea Identifyer" : "P6-255_ANT_2024_2025",
                                "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_antr2025_RINGS)

    # ANT 2024/25 - CHIRP2
    dict_UWB_antr2025_CHIRP2 = {"Campaign Folder"    : "2025_Antarctica_Polar6_rds_standard",
                                "Season Name"        : "ANT 2024/25",
                                "Radar System"       : "UWB",
                                "Platform"           : "Polar 6",
                                "Radar Specs"      : np.repeat("150-520 MHz bandwidth", 4),

                                "Segments"           : ["20241212_03",
                                                        "20241212_09",
                                                        "20241213_02",
                                                        "20241218_04"],

                                "Pangaea Event"     : ["P6-255_ANT_2024_2025_2412120801",
                                                        "P6-255_ANT_2024_2025_2412120801",
                                                        "P6-255_ANT_2024_2025_2412130901",
                                                        "P6-255_ANT_2024_2025_2412181401"],

                                "Campaign Name"      : "CHIRP 2",
                                "Campaign Name Long" : "Channel and Ice Rise Project in DML 2",
                                "Campaign PI"        : "Olaf Eisen (AWI), Reinhard Drews (University of Tübingen)",
                                "Grant ID"           : "AWI_PA_02159",
                                "Pangaea Identifyer" : "P6-255_ANT_2024_2025",
                                "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_antr2025_CHIRP2)

    # ANT 2024/25 - CHARISO
    dict_UWB_antr2025_CHARISO = {"Campaign Folder"   : "2025_Antarctica_Polar6_rds_standard",
                                "Season Name"        : "ANT 2024/25",
                                "Radar System"       : "UWB",
                                "Platform"           : "Polar 6",
                                "Radar Specs"      : np.repeat("150-520 MHz bandwidth", 9),

                                "Segments"           : ["20241214_01",
                                                        "20241214_02",
                                                        "20241214_03",
                                                        "20250101_04",
                                                        "20250102_01",
                                                        "20250102_03",
                                                        "20250105_01",
                                                        "20250105_03",
                                                        "20250107_02"],

                                "Pangaea Event"     : ["P6-255_ANT_2024_2025_2412141001",
                                                        "P6-255_ANT_2024_2025_2412141001",
                                                        "P6-255_ANT_2024_2025_2412141001",
                                                        "P6-255_ANT_2024_2025_2501012002",
                                                        "P6-255_ANT_2024_2025_2501022101",
                                                        "P6-255_ANT_2024_2025_2501022101",
                                                        "P6-255_ANT_2024_2025_2501052201",
                                                        "P6-255_ANT_2024_2025_2501052201",
                                                        "P6-255_ANT_2024_2025_2501072502"],

                                "Campaign Name"      : "CHARISO",
                                "Campaign Name Long" : "Isochrones in East Antarctic outflow regions (Swiss Antarctic Research Initiative)",
                                "Campaign PI"        : "Olaf Eisen (AWI), Johannes Sutter (University of Bern)",
                                "Grant ID"           : "AWI_PA_02145",
                                "Pangaea Identifyer" : "P6-255_ANT_2024_2025",
                                "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_antr2025_CHARISO)

    # ANT 2024/25 - SANAS
    dict_UWB_antr2025_SANAS = {"Campaign Folder"   : "2025_Antarctica_Polar6_rds_standard",
                                "Season Name"        : "ANT 2024/25",
                                "Radar System"       : "UWB",
                                "Platform"           : "Polar 6",
                                "Radar Specs"      : np.repeat("150-520 MHz bandwidth", 3),

                                "Segments"           : ["20250109_01",
                                                        "20250110_01",
                                                        "20250110_02"],

                                "Pangaea Event"     : ["P6-255_ANT_2024_2025_2501092601",
                                                        "P6-255_ANT_2024_2025_2501102701",
                                                        "P6-255_ANT_2024_2025_2501102802"],

                                "Campaign Name"      : "SANAS",
                                "Campaign Name Long" : "Structural Glaciological Analysis of Northern Antarctic Ice Shelf (SANAS) between Jelbart and Fimbul Shelves",
                                "Campaign PI"        : "Olaf Eisen (AWI), Sebastian Skatulla (University of Cape Town)",
                                "Grant ID"           : "AWI_PA_02161",
                                "Pangaea Identifyer" : "P6-255_ANT_2024_2025",
                                "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_antr2025_SANAS)

    # ANT 2024/25 - SISI
    dict_UWB_antr2025_SISI = {"Campaign Folder"      : "2025_Antarctica_Polar6_rds_standard",
                                "Season Name"        : "ANT 2024/25",
                                "Radar System"       : "UWB",
                                "Platform"           : "Polar 6",
                                "Radar Specs"      : np.repeat("150-520 MHz bandwidth", 1),
                                "Segments"           : ["20241222_03"],
                                "Pangaea Event"     : ["P6-255_ANT_2024_2025_2412221801"],
                                "Campaign Name"      : "SISI",
                                "Campaign Name Long" : "Detecting the Subsurface of blue Ice areas in Search for ancient Ice",
                                "Campaign PI"        : "Olaf Eisen (AWI), Francois Fripiat (Francois Fripiat)",
                                "Grant ID"           : "AWI_PA_02160",
                                "Pangaea Identifyer" : "P6-255_ANT_2024_2025",
                                "Radar Data DOI"     : "None"}

    list_UWB_ALL_dicts.append(dict_UWB_antr2025_SISI)


    ###################
    # All UWBM

    list_UWBM_ALL_dicts = []

    # ARK 2016 - Greenland
    dict_UWBM_arkr2016_Greenland    =  {"Campaign Folder"    : "2016_Greenland_Polar6_snow",
                                        "Season Name"        : "ARK 2016",
                                        "Radar System"       : "UWBM",
                                        "Platform"           : "Polar 6",
                                        "Radar Specs"        : np.repeat("2-18 GHz bandwidth", 28),

                                        "Segments"           : ["20160413_07",
                                                                "20160413_08",
                                                                "20160413_09",
                                                                "20160413_10",
                                                                "20160413_11",
                                                                "20160413_12",
                                                                "20160414_01",
                                                                "20160414_02",
                                                                "20160418_01",
                                                                "20160418_02",
                                                                "20160418_03",
                                                                "20160418_04",
                                                                "20160418_05",
                                                                "20160420_01",
                                                                "20160425_01",
                                                                "20160425_02",
                                                                "20160425_03",
                                                                "20160425_04",
                                                                "20160425_05",
                                                                "20160425_06",
                                                                "20160426_01",
                                                                "20160512_01",
                                                                "20160512_02",
                                                                "20160516_01",
                                                                "20160620_01",
                                                                "20160621_01",
                                                                "20160627_01",
                                                                "20160629_01"],

                                        "Pangaea Event"      : ["UWB_2016_1604130402",
                                                                "UWB_2016_1604130402",
                                                                "UWB_2016_1604130402",
                                                                "UWB_2016_1604130402",
                                                                "UWB_2016_1604130402",
                                                                "UWB_2016_1604130402",
                                                                "UWB_2016_1604140501",
                                                                "UWB_2016_1604140501",
                                                                "UWB_2016_1604180601",
                                                                "UWB_2016_1604180601",
                                                                "UWB_2016_1604180601",
                                                                "UWB_2016_1604180702",
                                                                "UWB_2016_1604180702",
                                                                "UWB_2016_1604200901",
                                                                "UWB_2016_1604251001",
                                                                "UWB_2016_1604251001",
                                                                "UWB_2016_1604251001",
                                                                "UWB_2016_1604251001",
                                                                "UWB_2016_1604251001",
                                                                "UWB_2016_1604251001",
                                                                "UWB_2016_1604261101",
                                                                "UWB_2016_1605121301",
                                                                "UWB_2016_1605121301",
                                                                "UWB_2016_1605161401",
                                                                "UWB_2016_1606201601",
                                                                "UWB_2016_1606211701",
                                                                "UWB_2016_1606271801",
                                                                "UWB_2016_1606291901"],

                                        "Campaign Name"      : "UWB/M Testing",
                                        "Campaign Name Long" : "Testing of the new UWB and UWB/M radars",
                                        "Campaign PI"        : "Heinrich Miller (AWI)",
                                        "Grant ID"           : "AWI_PA_02078",
                                        "Pangaea Identifyer" : "P6-201_UWB_2016",
                                        "Radar Data DOI"     : "None"}

    list_UWBM_ALL_dicts.append(dict_UWBM_arkr2016_Greenland)

    # ARK 2017 - PAMARCMIP
    dict_UWBM_arkr2017_PAMARCMIP    =  {"Campaign Folder"    : "2017_Arctic_PAMARCMIP_Polar5_snow",
                                        "Season Name"        : "ARK 2017",
                                        "Radar System"       : "UWBM",
                                        "Platform"           : "Polar 5",
                                        "Radar Specs"      : np.repeat("2-18 GHz bandwidth", 16),

                                        "Segments"           : ["20170303_01",
                                                                "20170303_02",
                                                                "20170321_01",
                                                                "20170324_01",
                                                                "20170324_02",
                                                                "20170328_01",
                                                                "20170328_02",
                                                                "20170328_03",
                                                                "20170330_01",
                                                                "20170330_03",
                                                                "20170402_01",
                                                                "20170404_01",
                                                                "20170404_02",
                                                                "20170406_01",
                                                                "20170408_01",
                                                                "20170410_01"],

                                        "Pangaea Event"      : ["P5_205_PAMARCMIP_2017_1703030101",
                                                                "P5_205_PAMARCMIP_2017_1703030101",
                                                                "P5_205_PAMARCMIP_2017_1703210501",
                                                                "P5_205_PAMARCMIP_2017_1703240801",
                                                                "P5_205_PAMARCMIP_2017_1703240801",
                                                                "P5_205_PAMARCMIP_2017_1703280901",
                                                                "P5_205_PAMARCMIP_2017_1703280901",
                                                                "P5_205_PAMARCMIP_2017_1703280901",
                                                                "P5_205_PAMARCMIP_2017_1703301101",
                                                                "P5_205_PAMARCMIP_2017_1703311201",
                                                                "P5_205_PAMARCMIP_2017_1704021401",
                                                                "P5_205_PAMARCMIP_2017_1704041501",
                                                                "P5_205_PAMARCMIP_2017_1704041501",
                                                                "P5_205_PAMARCMIP_2017_1704061801",
                                                                "P5_205_PAMARCMIP_2017_1704082001",
                                                                "P5_205_PAMARCMIP_2017_1704102101"],

                                        "Campaign Name"      : "PAMARCMIP 2017",
                                        "Campaign Name Long" : "Polar Airborne Measurements and Arctic regional Climate",
                                        "Campaign PI"        : "Andreas Herber (AWI)",
                                        "Grant ID"           : "AWI_PA_02086",
                                        "Pangaea Identifyer" : "P5-205_PAMARCMIP_2017",
                                        "Radar Data DOI"     : "None"}

    list_UWBM_ALL_dicts.append(dict_UWBM_arkr2017_PAMARCMIP)


    # ARK 2018 - RESURV79
    dict_UWBM_antr2018_RESURV79     =  {"Campaign Folder"    : "2018_Greenland_SURV79_Polar6_snow",
                                        "Season Name"        : "ARK 2018",
                                        "Radar System"       : "UWBM",
                                        "Platform"           : "Polar 6",
                                        "Radar Specs"        : np.repeat("2-18 GHz bandwidth", 24),

                                        "Segments"           : ["20180329_01",
                                                                "20180329_02",
                                                                "20180404_01",
                                                                "20180410_01",
                                                                "20180411_01",
                                                                "20180413_01",
                                                                "20180414_01",
                                                                "20180415_02",
                                                                "20180415_03",
                                                                "20180415_04",
                                                                "20180417_01",
                                                                "20180418_01",
                                                                "20180418_02",
                                                                "20180418_03",
                                                                "20180420_01",
                                                                "20180422_01",
                                                                "20180422_02",
                                                                "20180423_01",
                                                                "20180423_02",
                                                                "20180423_03",
                                                                "20180423_04",
                                                                "20180423_05",
                                                                "20180423_06",
                                                                "20180423_07"],

                                        "Pangaea Event"     :  ["P6_211_RESURV79_2018_1803290101",
                                                                "P6_211_RESURV79_2018_1803290101",
                                                                "P6_211_RESURV79_2018_1804040201",
                                                                "P6_211_RESURV79_2018_1804100301",
                                                                "P6_211_RESURV79_2018_1804110401",
                                                                "P6_211_RESURV79_2018_1804130501",
                                                                "P6_211_RESURV79_2018_1804140601",
                                                                "P6_211_RESURV79_2018_1804150701",
                                                                "P6_211_RESURV79_2018_1804150701",
                                                                "P6_211_RESURV79_2018_1804150701",
                                                                "P6_211_RESURV79_2018_1804170801",
                                                                "P6_211_RESURV79_2018_1804180901",
                                                                "P6_211_RESURV79_2018_1804180901",
                                                                "P6_211_RESURV79_2018_1804180901",
                                                                "P6_211_RESURV79_2018_1804201001",
                                                                "P6_211_RESURV79_2018_1804221101",
                                                                "P6_211_RESURV79_2018_1804221101",
                                                                "P6_211_RESURV79_2018_1804231201",
                                                                "P6_211_RESURV79_2018_1804231201",
                                                                "P6_211_RESURV79_2018_1804231201",
                                                                "P6_211_RESURV79_2018_1804231201",
                                                                "P6_211_RESURV79_2018_1804231201",
                                                                "P6_211_RESURV79_2018_1804231201",
                                                                "P6_211_RESURV79_2018_1804231201"],

                                        "Campaign Name"      : "RESURV79",
                                        "Campaign Name Long" : "Re-survey of the 79° Glacier",
                                        "Campaign PI"        : "Angelika Humbert (AWI)",
                                        "Grant ID"           : "AWI_PA_02095",
                                        "Pangaea Identifyer" : "P6-211_RESURV79_2018",
                                        "Radar Data DOI"     : "None"}

    list_UWBM_ALL_dicts.append(dict_UWBM_antr2018_RESURV79)


    # ARK 2019 - ICEBIRD
    dict_UWBM_arkr2019_ICEBIRD      =  {"Campaign Folder"    : "2019_Arctic_ICEBIRD_Polar6_snow",
                                        "Season Name"        : "ARK 2019",
                                        "Radar System"       : "UWBM",
                                        "Platform"           : "Polar 6",
                                        "Radar Specs"        : np.repeat("2-18 GHz bandwidth", 12),

                                        "Segments"           : ["20190401_01",
                                                                "20190402_01",
                                                                "20190402_02",
                                                                "20190405_02",
                                                                "20190407_01",
                                                                "20190407_02",
                                                                "20190407_04",
                                                                "20190408_01",
                                                                "20190408_02",
                                                                "20190410_01",
                                                                "20190410_02",
                                                                "20190410_03"],

                                        "Pangaea Event"     : ["P6_217_ICEBIRD_2019_1904010801",
                                                                "P6_217_ICEBIRD_2019_1904020901",
                                                                "P6_217_ICEBIRD_2019_1904020901",
                                                                "P6_217_ICEBIRD_2019_1904051001",
                                                                "P6_217_ICEBIRD_2019_1904071201",
                                                                "P6_217_ICEBIRD_2019_1904071201",
                                                                "P6_217_ICEBIRD_2019_1904071201",
                                                                "P6_217_ICEBIRD_2019_1904081301",
                                                                "P6_217_ICEBIRD_2019_1904081301",
                                                                "P6_217_ICEBIRD_2019_1904101401",
                                                                "P6_217_ICEBIRD_2019_1904101401",
                                                                "P6_217_ICEBIRD_2019_1904101401"],

                                        "Campaign Name"      : "IceBird 2019 Winter",
                                        "Campaign Name Long" : "ICESat-2 Validation Data Acquisition",
                                        "Campaign PI"        : "Stefan Hendricks (AWI)",
                                        "Grant ID"           : "AWI_PA_02108",
                                        "Pangaea Identifyer" : "P6-217_ICEBIRD_2019",
                                        "Radar Data DOI"     : "None"}

    list_UWBM_ALL_dicts.append(dict_UWBM_arkr2019_ICEBIRD)


    # ANT 2022/23 - ANTSI (Atka Bay)
    dict_UWBM_antr2023_ANTSI_AtkaBay = {"Campaign Folder"    : "2022_Antarctica_Polar5_snow",
                                        "Season Name"        : "ANT 2022/23",
                                        "Radar System"       : "UWBM",
                                        "Platform"           : "Polar 5",
                                        "Radar Specs"        : np.repeat("2-18 GHz bandwidth", 38),

                                        "Segments"           : ["20221205_01",
                                                                "20221205_02",
                                                                "20221205_03",
                                                                "20221205_04",
                                                                "20221205_05",
                                                                "20221205_06",
                                                                "20221205_07",
                                                                "20221205_08",
                                                                "20221205_09",
                                                                "20221205_10",
                                                                "20221205_11",
                                                                "20221205_12",
                                                                "20221205_13",
                                                                "20221205_14",
                                                                "20221205_15",
                                                                "20221205_16",
                                                                "20221205_17",
                                                                "20221205_18",
                                                                "20221205_19",
                                                                "20221205_20",
                                                                "20221205_21",
                                                                "20221205_22",
                                                                "20221205_23",
                                                                "20221205_24",
                                                                "20221205_25",
                                                                "20221205_26",
                                                                "20221205_27",
                                                                "20221205_28",
                                                                "20221205_29",
                                                                "20221205_30",
                                                                "20221205_31",
                                                                "20221205_32",
                                                                "20221205_33",
                                                                "20221205_34",
                                                                "20221205_35",
                                                                "20221205_36",
                                                                "20221205_37",
                                                                "20221205_38"],

                                        "Pangaea Event"     : np.repeat("P5_238_ANTSI_2022_2212051101", 38),

                                        "Campaign Name"      : "ANTSI (Atka Bay)",
                                        "Campaign Name Long" : "Antarctic Sea Ice: Thickness, Melt Ponding, and Ice Shelf Interaction",
                                        "Campaign PI"        : "Christian Haas (AWI)",
                                        "Grant ID"           : "AWI_PA_02135",
                                        "Pangaea Identifyer" : "P5-238_ANTSI_2022",
                                        "Radar Data DOI"     : "https://doi.org/10.1594/PANGAEA.981151"}

    list_UWBM_ALL_dicts.append(dict_UWBM_antr2023_ANTSI_AtkaBay)

    # ANT 2022/23 - ANTSI (Blue Ice - Troll Station)
    dict_UWBM_antr2023_ANTSI_Troll   = {"Campaign Folder"    : "2022_Antarctica_Polar5_snow",
                                        "Season Name"        : "ANT 2022/23",
                                        "Radar System"       : "UWBM",
                                        "Platform"           : "Polar 5",
                                        "Radar Specs"        : np.repeat("2-18 GHz bandwidth", 1),
                                        "Segments"           : ["20221203_02"],
                                        "Pangaea Event"      : ["P5_238_ANTSI_2022_2212030901"],
                                        "Campaign Name"      : "ANTSI (Blue Ice - Troll Station)",
                                        "Campaign Name Long" : "Antarctic Sea Ice: Thickness, Melt Ponding, and Ice Shelf Interaction",
                                        "Campaign PI"        : "Christian Haas (AWI)",
                                        "Grant ID"           : "AWI_PA_02135",
                                        "Pangaea Identifyer" : "P5-238_ANTSI_2022",
                                        "Radar Data DOI"     : "None"}

    list_UWBM_ALL_dicts.append(dict_UWBM_antr2023_ANTSI_Troll)

    # ANT 2022/23 - Kottas-Aero
    dict_UWBM_antr2023_Kottas_Aero   = {"Campaign Folder"    : "2022_Antarctica_Polar5_snow",
                                        "Season Name"        : "ANT 2022/23",
                                        "Radar System"       : "UWBM",
                                        "Platform"           : "Polar 5",
                                        "Radar Specs"        : np.repeat("2-18 GHz bandwidth", 3),
                                        "Segments"           : ["20221204_01",                   "20221205_99",                 "20221209_01"],
                                        "Pangaea Event"      : ["P5_238_ANTSI_2022_2212041001", "P5_238_ANTSI_2022_2212051202", "P5_238_ANTSI_2022_2212091301"],
                                        "Campaign Name"      : "Kottas-Aero",
                                        "Campaign Name Long" : "Distributed surface mass balance and elevation along the Kottaspegel traverse",
                                        "Campaign PI"        : "Olaf Eisen (AWI)",
                                        "Grant ID"           : "AWI_PA_02137",
                                        "Pangaea Identifyer" : "P5-238_ANTSI_2022",
                                        "Radar Data DOI"     : "None"}

    list_UWBM_ALL_dicts.append(dict_UWBM_antr2023_Kottas_Aero)

    # ARK 2023 - IceBird Winter
    dict_UWBM_arkr2023_IceBird       =  {"Campaign Folder"    : "2023_Arctic_Polar6",
                                        "Season Name"        : "ARK 2023",
                                        "Radar System"       : "UWBM",
                                        "Platform"           : "Polar 6",
                                        "Radar Specs"        : np.repeat("2-18 GHz bandwidth", 8),

                                        "Segments"           : ["20230327_01",
                                                                "20230331_01",
                                                                "20230402_01",
                                                                "20230404_01",
                                                                "20230405_01",
                                                                "20230408_01",
                                                                "20230419_01",
                                                                "20230420_01"],

                                        "Pangaea Event"      : ["P6_239_IceBird_Winter_2023_2303270201",
                                                                "P6_239_IceBird_Winter_2023_2303310301",
                                                                "P6_239_IceBird_Winter_2023_2304020401",
                                                                "P6_239_IceBird_Winter_2023_2304040501",
                                                                "P6_239_IceBird_Winter_2023_2304050601",
                                                                "P6_239_IceBird_Winter_2023_2304080701",
                                                                "P6_239_IceBird_Winter_2023_2304190801",
                                                                "P6_239_IceBird_Winter_2023_2304200901"],

                                        "Campaign Name"      : "IceBird-W23",
                                        "Campaign Name Long" : "IceBird Winter 2023 sea ice observation program",
                                        "Campaign PI"        : "Christian Haas (AWI)",
                                        "Grant ID"           : "AWI_PA_02139",
                                        "Pangaea Identifyer" : "P6-239_IceBird_Winter_2023",
                                        "Radar Data DOI"     : "None"}

    list_UWBM_ALL_dicts.append(dict_UWBM_arkr2023_IceBird)

    # ARK 2024 - IceBird Winter
    dict_UWBM_arkr2024_IceBird      =  {"Campaign Folder"    : "2024_Arctic_ICEBIRD_Polar5_snow",
                                        "Season Name"        : "ARK 2024",
                                        "Radar System"       : "UWBM",
                                        "Platform"           : "Polar 5",
                                        "Radar Specs"      : np.repeat("2-18 GHz bandwidth", 36),

                                        "Segments"           : ["20240406_01",
                                                                "20240406_02",
                                                                "20240406_03",
                                                                "20240411_01",
                                                                "20240412_01",
                                                                "20240415_01",
                                                                "20240415_02",
                                                                "20240415_03",
                                                                "20240415_04",
                                                                "20240415_05",
                                                                "20240416_01",
                                                                "20240416_02",
                                                                "20240416_03",
                                                                "20240416_04",
                                                                "20240418_01",
                                                                "20240418_02",
                                                                "20240418_03",
                                                                "20240418_04",
                                                                "20240418_05",
                                                                "20240420_01",
                                                                "20240420_02",
                                                                "20240420_03",
                                                                "20240421_01",
                                                                "20240421_02",
                                                                "20240421_03",
                                                                "20240421_04",
                                                                "20240426_01",
                                                                "20240426_02",
                                                                "20240427_01",
                                                                "20240429_01",
                                                                "20240429_02",
                                                                "20240502_01",
                                                                "20240502_02",
                                                                "20240502_03",
                                                                "20240502_04",
                                                                "20240506_01"],

                                        "Pangaea Event"     : ["P5-248_IceBird_Winter_2024_2404060101",
                                                                "P5-248_IceBird_Winter_2024_2404060101",
                                                                "P5-248_IceBird_Winter_2024_2404060101",
                                                                "P5-248_IceBird_Winter_2024_2404110201",
                                                                "P5-248_IceBird_Winter_2024_2404120301",
                                                                "P5-248_IceBird_Winter_2024_2404150401",
                                                                "P5-248_IceBird_Winter_2024_2404150401",
                                                                "P5-248_IceBird_Winter_2024_2404150401",
                                                                "P5-248_IceBird_Winter_2024_2404150401",
                                                                "P5-248_IceBird_Winter_2024_2404150401",
                                                                "P5-248_IceBird_Winter_2024_2404160501",
                                                                "P5-248_IceBird_Winter_2024_2404160501",
                                                                "P5-248_IceBird_Winter_2024_2404160501",
                                                                "P5-248_IceBird_Winter_2024_2404160501",
                                                                "P5-248_IceBird_Winter_2024_2404180601",
                                                                "P5-248_IceBird_Winter_2024_2404180601",
                                                                "P5-248_IceBird_Winter_2024_2404180601",
                                                                "P5-248_IceBird_Winter_2024_2404180601",
                                                                "P5-248_IceBird_Winter_2024_2404180601",
                                                                "P5-248_IceBird_Winter_2024_2404200701",
                                                                "P5-248_IceBird_Winter_2024_2404200802",
                                                                "P5-248_IceBird_Winter_2024_2404200802",
                                                                "P5-248_IceBird_Winter_2024_2404210901",
                                                                "P5-248_IceBird_Winter_2024_2404210901",
                                                                "P5-248_IceBird_Winter_2024_2404210901",
                                                                "P5-248_IceBird_Winter_2024_2404210901",
                                                                "P5-248_IceBird_Winter_2024_2404261001",
                                                                "P5-248_IceBird_Winter_2024_2404261001",
                                                                "P5-248_IceBird_Winter_2024_2404271101",
                                                                "P5-248_IceBird_Winter_2024_2404291201",
                                                                "P5-248_IceBird_Winter_2024_2404291201",
                                                                "P5-248_IceBird_Winter_2024_2405021301",
                                                                "P5-248_IceBird_Winter_2024_2405021301",
                                                                "P5-248_IceBird_Winter_2024_2405021301",
                                                                "P5-248_IceBird_Winter_2024_2405021301",
                                                                "P5-248_IceBird_Winter_2024_2405061401"],

                                        "Campaign Name"      : "IceBird-CAN24",
                                        "Campaign Name Long" : "IceBird sea ice thickness surveys in the Canadian Arctic in 2024",
                                        "Campaign PI"        : "Christian Haas (AWI)",
                                        "Grant ID"           : "AWI_PA_02150",
                                        "Pangaea Identifyer" : "P5-248_IceBird_Winter_2024",
                                        "Radar Data DOI"     : "None"}

    list_UWBM_ALL_dicts.append(dict_UWBM_arkr2024_IceBird)

    # ARK 2025 - ICEBIRD Winter
    dict_UWBM_arkr2025_ICEBIRD      =  {"Campaign Folder"    : "2025_Arctic_ICEBIRD_Polar6_snow",
                                        "Season Name"        : "ARK 2025",
                                        "Radar System"       : "UWBM",
                                        "Platform"           : "Polar 6",
                                        "Radar Specs"        : np.repeat("2-18 GHz bandwidth", 78),

                                        "Segments"           : ["20250313_01",
                                                                "20250320_01",
                                                                "20250322_01",
                                                                "20250322_02",
                                                                "20250323_01",
                                                                "20250323_02",
                                                                "20250325_01",
                                                                "20250329_01",
                                                                "20250329_02",
                                                                "20250329_03",
                                                                "20250329_04",
                                                                "20250329_05",
                                                                "20250329_06",
                                                                "20250329_07",
                                                                "20250329_08",
                                                                "20250331_01",
                                                                "20250331_02",
                                                                "20250331_03",
                                                                "20250331_04",
                                                                "20250331_05",
                                                                "20250331_06",
                                                                "20250331_07",
                                                                "20250331_08",
                                                                "20250331_09",
                                                                "20250402_01",
                                                                "20250402_02",
                                                                "20250402_03",
                                                                "20250402_04",
                                                                "20250402_05",
                                                                "20250403_01",
                                                                "20250403_02",
                                                                "20250403_03",
                                                                "20250403_04",
                                                                "20250403_05",
                                                                "20250403_06",
                                                                "20250403_07",
                                                                "20250403_08",
                                                                "20250403_09",
                                                                "20250403_10",
                                                                "20250403_11",
                                                                "20250405_01",
                                                                "20250405_02",
                                                                "20250405_03",
                                                                "20250405_04",
                                                                "20250406_01",
                                                                "20250406_02",
                                                                "20250406_03",
                                                                "20250406_04",
                                                                "20250406_05",
                                                                "20250406_06",
                                                                "20250406_07",
                                                                "20250407_02",
                                                                "20250407_03",
                                                                "20250407_04",
                                                                "20250407_05",
                                                                "20250407_06",
                                                                "20250407_07",
                                                                "20250407_08",
                                                                "20250407_09",
                                                                "20250407_10",
                                                                "20250413_02",
                                                                "20250414_01",
                                                                "20250414_02",
                                                                "20250414_03",
                                                                "20250414_04",
                                                                "20250414_05",
                                                                "20250414_06",
                                                                "20250417_01",
                                                                "20250417_02",
                                                                "20250417_03",
                                                                "20250418_01",
                                                                "20250419_03",
                                                                "20250420_01",
                                                                "20250421_02",
                                                                "20250421_03",
                                                                "20250421_04",
                                                                "20250421_05",
                                                                "20250421_06"],    

                                        "Pangaea Event"      : ["P6-257_IceBird_Winter_2025_2503130101",
                                                                "P6-257_IceBird_Winter_2025_2503200301",
                                                                "P6-257_IceBird_Winter_2025_2503220401",
                                                                "P6-257_IceBird_Winter_2025_2503220401",
                                                                "P6-257_IceBird_Winter_2025_2503230501",
                                                                "P6-257_IceBird_Winter_2025_2503230501",
                                                                "P6-257_IceBird_Winter_2025_2503250601",
                                                                "P6-257_IceBird_Winter_2025_2503290801",
                                                                "P6-257_IceBird_Winter_2025_2503290801",
                                                                "P6-257_IceBird_Winter_2025_2503290801",
                                                                "P6-257_IceBird_Winter_2025_2503290801",
                                                                "P6-257_IceBird_Winter_2025_2503290801",
                                                                "P6-257_IceBird_Winter_2025_2503290801",
                                                                "P6-257_IceBird_Winter_2025_2503290801",
                                                                "P6-257_IceBird_Winter_2025_2503290801",
                                                                "P6-257_IceBird_Winter_2025_2503311001",
                                                                "P6-257_IceBird_Winter_2025_2503311001",
                                                                "P6-257_IceBird_Winter_2025_2503311001",
                                                                "P6-257_IceBird_Winter_2025_2503311001",
                                                                "P6-257_IceBird_Winter_2025_2503311001",
                                                                "P6-257_IceBird_Winter_2025_2503311001",
                                                                "P6-257_IceBird_Winter_2025_2503311001",
                                                                "P6-257_IceBird_Winter_2025_2503311001",
                                                                "P6-257_IceBird_Winter_2025_2503311001",
                                                                "P6-257_IceBird_Winter_2025_2504021101",
                                                                "P6-257_IceBird_Winter_2025_2504021101",
                                                                "P6-257_IceBird_Winter_2025_2504021101",
                                                                "P6-257_IceBird_Winter_2025_2504021101",
                                                                "P6-257_IceBird_Winter_2025_2504021101",
                                                                "P6-257_IceBird_Winter_2025_2504031201",
                                                                "P6-257_IceBird_Winter_2025_2504031201",
                                                                "P6-257_IceBird_Winter_2025_2504031201",
                                                                "P6-257_IceBird_Winter_2025_2504031201",
                                                                "P6-257_IceBird_Winter_2025_2504031201",
                                                                "P6-257_IceBird_Winter_2025_2504031201",
                                                                "P6-257_IceBird_Winter_2025_2504031201",
                                                                "P6-257_IceBird_Winter_2025_2504031201",
                                                                "P6-257_IceBird_Winter_2025_2504031201",
                                                                "P6-257_IceBird_Winter_2025_2504031201",
                                                                "P6-257_IceBird_Winter_2025_2504031201",
                                                                "P6-257_IceBird_Winter_2025_2504051301",
                                                                "P6-257_IceBird_Winter_2025_2504051301",
                                                                "P6-257_IceBird_Winter_2025_2504051301",
                                                                "P6-257_IceBird_Winter_2025_2504051301",
                                                                "P6-257_IceBird_Winter_2025_2504061401",
                                                                "P6-257_IceBird_Winter_2025_2504061401",
                                                                "P6-257_IceBird_Winter_2025_2504061401",
                                                                "P6-257_IceBird_Winter_2025_2504061401",
                                                                "P6-257_IceBird_Winter_2025_2504061401",
                                                                "P6-257_IceBird_Winter_2025_2504061401",
                                                                "P6-257_IceBird_Winter_2025_2504061401",
                                                                "P6-257_IceBird_Winter_2025_2504071601",
                                                                "P6-257_IceBird_Winter_2025_2504071601",
                                                                "P6-257_IceBird_Winter_2025_2504071601",
                                                                "P6-257_IceBird_Winter_2025_2504071601",
                                                                "P6-257_IceBird_Winter_2025_2504071601",
                                                                "P6-257_IceBird_Winter_2025_2504071601",
                                                                "P6-257_IceBird_Winter_2025_2504071601",
                                                                "P6-257_IceBird_Winter_2025_2504071601",
                                                                "P6-257_IceBird_Winter_2025_2504071601",
                                                                "P6-257_IceBird_Winter_2025_2504131701",
                                                                "P6-257_IceBird_Winter_2025_2504141801",
                                                                "P6-257_IceBird_Winter_2025_2504141801",
                                                                "P6-257_IceBird_Winter_2025_2504141801",
                                                                "P6-257_IceBird_Winter_2025_2504141801",
                                                                "P6-257_IceBird_Winter_2025_2504141801",
                                                                "P6-257_IceBird_Winter_2025_2504141801",
                                                                "P6-257_IceBird_Winter_2025_2504171901",
                                                                "P6-257_IceBird_Winter_2025_2504171901",
                                                                "P6-257_IceBird_Winter_2025_2504171901",
                                                                "P6-257_IceBird_Winter_2025_2504182001",
                                                                "P6-257_IceBird_Winter_2025_2504192101",
                                                                "P6-257_IceBird_Winter_2025_2504202201",
                                                                "P6-257_IceBird_Winter_2025_2504212301",
                                                                "P6-257_IceBird_Winter_2025_2504212301",
                                                                "P6-257_IceBird_Winter_2025_2504212301",
                                                                "P6-257_IceBird_Winter_2025_2504212301",
                                                                "P6-257_IceBird_Winter_2025_2504212301"],

                                        "Campaign Name"      : "IceBird-CAN25",
                                        "Campaign Name Long" : "IceBird airborne snow and sea ice surveys in the Canadian Arctic (in winter 2025)",
                                        "Campaign PI"        : "Christian Haas (AWI), Luisa Wagner (AWI), Richard Kelly (University of Waterloo, ON, Canada)",
                                        "Grant ID"           : "AWI_PA_02162",
                                        "Pangaea Identifyer" : "P6-257_IceBird_Winter_2025",
                                        "Radar Data DOI"     : "None"}

    list_UWBM_ALL_dicts.append(dict_UWBM_arkr2025_ICEBIRD)

    return dict_ASIRAS, dict_SNOW, dict_ACCU, dict_EMR_ANT, dict_EMR_ARK, list_UWB_ALL_dicts, list_UWBM_ALL_dicts