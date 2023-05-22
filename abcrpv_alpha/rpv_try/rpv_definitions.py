# Sparticles
SPARTICLES=["B","W^+", "W^0","G","H^+","H^0","q","d","u","l","nu","e","b_L","t_L","t","b","tau_L","nu_tau","tau"]
# Notation used in for sparticles in code:
# "B",      = Bino
# "W^+",    = Charged Wino     
# "W^0",    = Neutral Wino      
# "G",      = Gluino 
# "H^+",    = Charged Higgsino     
# "H^0",    = Neutral Higgsino     
# "q",      = u_L, d_L, c_L, s_L
# "d",      = d_R, s_R
# "u",      = u_R, c_R
# "l",      = e_L, mu_L 
# "nu",     = nu_e, nu_mu 
# "e",      = e_R, mu_R
# "t_L",    = t_L     
# "b_L",    = b_L     
# "t",      = t_R 
# "b",      = b_R 
# "tau_L",  = tau_L
# "tau"     = tau_R
# "nu_tau", = nu_tau     

# Final State Object
FINAL_STATE=["v","J","3","t","b","j","L","T","l","E"]
# Notation used in for fiinal state object in code (one-char syntax):
# "v" = Massive Bosons  (W,Z,H)  
# "J" = jet             (u, d, c, s, t, b )         
# "3" = 3rd gen jet     (t, b )         
# "t" = top             (t)         
# "b" = bottom jet      (b)                 
# "j" = light jet       (u, d, c, s)                 
# "L" = charged leptons (e, mu, tau )                     
# "l" = light lepton    (e, mu )                     
# "T" = tau             (tau)         
# "E" = MET             (nu)         

# sparticles that are relevant to RPV decays
RPV_sparticles = ["q","d","u","l","nu","e","b_L","t_L","t","b","tau_L","nu_tau","tau"]

#Categories for each RPV couplings 
LLE_CAT = ["L L E", "L_3 L E", "L L E_3", "L_3 L E_3"]
LQD_CAT = ["L Q D", "L Q D_3", "L Q_3 D", "L Q_3 D_3", "L_3 Q D", "L_3 Q D_3", "L_3 Q_3 D", "L_3 Q_3 D_3"]
UDD_CAT = ["U D D","U D_3 D","U D_3 D_3","U_3 D D","U_3 D_3 D","U_3 D_3 D_3"]

#"vertex" according to the Categories
LLE_STATES= (
[ "l"      ,   "nu"         ,   "e"    ,   "L L E"      ], 
[ "l"      ,   "nu_tau"     ,   "e"    ,   "L_3 L E"    ], 
[ "tau_L"  ,   "nu"         ,   "e"    ,   "L_3 L E"    ], 
[ "l"      ,   "nu"         ,   "tau"  ,   "L L E_3"    ], 
[ "l"      ,   "nu_tau"     ,   "tau"  ,   "L_3 L E_3"  ], 
[ "tau_L"  ,   "nu"         ,   "tau"  ,   "L_3 L E_3"  ])

LQD_STATES=(
["l"       , "q"     , "d"   , "L Q D"        ]  ,  ["nu"      , "q"    , "d"   , "L Q D"          ],
["l"       , "q"     , "b"   , "L Q D_3"      ]  ,  ["nu"      , "q"    , "b"   , "L Q D_3"        ],
["l"       , "t_L"   , "d"   , "L Q_3 D"      ]  ,  ["nu"      , "b_L"  , "d"   , "L Q_3 D"        ],
["l"       , "t_L"   , "b"   , "L Q_3 D_3"    ]  ,  ["nu"      , "b_L"  , "b"   , "L Q_3 D_3"      ],
["tau_L"   , "q"     , "d"   , "L_3 Q D"      ]  ,  ["nu_tau"  , "q"    , "d"   , "L_3 Q D"        ],
["tau_L"   , "q"     , "b"   , "L_3 Q D_3"    ]  ,  ["nu_tau"  , "q"    , "b"   , "L_3 Q D_3"      ],
["tau_L"   , "t_L"   , "d"   , "L_3 Q_3 D"    ]  ,  ["nu_tau"  , "b_L"  , "d"   , "L_3 Q_3 D"      ],
["tau_L"   , "t_L"   , "b"   , "L_3 Q_3 D_3"  ]  ,  ["nu_tau"  , "b_L"  , "b"   , "L_3 Q_3 D_3"    ])

UDD_STATES=(
["u"       , "d"    , "d"   ,   "U D D"        ],
["u"       , "d"    , "b"   ,   "U D_3 D"      ],
["u"       , "b"    , "b"   ,   "U D_3 D_3"    ],
["t"       , "d"    , "d"   ,   "U_3 D D"      ],
["t"       , "d"    , "b"   ,   "U_3 D_3 D"    ],
["t"       , "b"    , "b"   ,   "U_3 D_3 D_3"  ])

sig_order = ["v","J","3","t","b","j","L","T","l","E"]
