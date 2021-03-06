#! /usr/bin/env python                                                                                                                           
# -*- coding: utf-8 -*-                                                                                                                          
 
def lit_data():
    class fluxes():
        def __init__(self):
            self.totalFluxes = {}
            self.title = []
    liq10 = fluxes()
    liq10.totalFluxes = {
             'Fram' : {'Inflow' : 3.9, 'Outflow' : -5.5, 'Total flow' : -1.6 , 'FWFlux' : -47.6, 'HFlux':24.4} ,
             'Fram1' : {'Inflow' : 3.9, 'Outflow' : -5.5, 'Total flow' : -1.6, 'FWFlux': -47.6, 'HFlux':24.4} ,
             'Fram2' : {'Inflow' : 3.9, 'Outflow' : -5.5, 'Total flow' : -1.6, 'FWFlux' : -47.6, 'HFlux':24.4} ,
             'Barents' : {'Inflow' : -1.2, 'Outflow' : 4.1, 'Total flow' :2.9, 'FWFlux':-8.1, 'HFlux':71.6} ,
             'Barents1' : {'Inflow' : -1.2, 'Outflow' : 4.1, 'Total flow' :2.9, 'FWFlux':-8.1, 'HFlux':71.6} ,
             'Barents2' : {'Inflow' : -1.2, 'Outflow' : 4.1, 'Total flow' :2.9, 'FWFlux' :-8.1, 'HFlux':71.6} ,
             'Bering': {'Inflow' : 1.3, 'Outflow' : 0.0, 'Total flow' :1.3, 'FWFlux' :95.2,'HFlux':5.4} ,
             'Davis': {'Inflow' : 0.6, 'Outflow' : -3.2, 'Total flow' :-2.6, 'FWFlux' :-123.8, 'HFlux':12.6},
             'Davis2': {'Inflow' : 0.6, 'Outflow' : -3.2, 'Total flow' :-2.6, 'FWFlux' :-123.8, 'HFlux':12.6},
             'Davis1': {'Inflow' : 0.6, 'Outflow' : -3.2, 'Total flow' :-2.6, 'FWFlux' :-123.8, 'HFlux':12.6}
            }
    liq10.color = 'blue'
    liq10.title = "Lique10"
    
    core2 = fluxes()
    core2.totalFluxes =  {
             'Fram' : {'Inflow' : 0, 'Outflow' : 0, 'Total flow' : -1.95, 'FWFlux' : 0.0, 'HFlux':53.96,'Error':2.3} ,
             'Fram1' : {'Inflow' : 0, 'Outflow' : 0, 'Total flow' : -1.95, 'FWFlux' : 0.0, 'HFlux':53.96,'Error':2.3} ,
             'Fram2' : {'Inflow' : 0, 'Outflow' : 0, 'Total flow' : -1.95, 'FWFlux' : 0.0, 'HFlux':53.96,'Error':1.12} ,
             'Barents' : {'Inflow' : 0, 'Outflow' : 0, 'Total flow' : 2.53, 'FWFlux' : 0.0, 'HFlux':22.98,'Error':1.06} ,
             'Barents1' : {'Inflow' : 0, 'Outflow' : 0, 'Total flow' : 2.53, 'FWFlux' : 0.0, 'HFlux':22.98,'Error':1.06} ,
             'Barents2' : {'Inflow' : 0, 'Outflow' : 0, 'Total flow' : 2.53, 'FWFlux' : 0.0, 'HFlux':22.98,'Error':1.06} ,
             'Bering': {'Inflow' : 0, 'Outflow' :0, 'Total flow' : 0.99, 'FWFlux' : 0.0, 'HFlux': 3.39,'Error':0.3} ,
             'Davis': {'Inflow' : 0, 'Outflow' : 0, 'Total flow' :-1.75, 'FWFlux' : 0.0, 'HFlux':13.44,'Error':1.2},
             'Davis2': {'Inflow' : 0, 'Outflow' : 0, 'Total flow' :-1.75, 'FWFlux' : 0.0, 'HFlux':13.44,'Error':1.2},
             'Davis1': {'Inflow' : 0, 'Outflow' : 0, 'Total flow' :-1.75, 'FWFlux' : 0.0, 'HFlux':13.44,'Error':1.2} 
            }
    core2.title = "CORE2"
    core2.color = 'orange'

    obs = fluxes()
    obs.totalFluxes =   {'Fram' : {'Inflow' : 0, 'Outflow' : 0, 'Total flow' : -2.0, 'FWFlux' : 0.0, 'HFlux':70.0,'Error':-2.7} ,
             'Fram1' : {'Inflow' : 0, 'Outflow' : 0, 'Total flow' : -2.0, 'FWFlux' : 0.0, 'HFlux':70.0,'Error':-2.7} ,
             'Fram2' : {'Inflow' : 0, 'Outflow' : 0, 'Total flow' : -2.0, 'FWFlux' : 0.0, 'HFlux':70.0,'Error':-2.7} ,
             'Barents' : {'Inflow' : -1.2, 'Outflow' : 4.1, 'Total flow' : 2.15 , 'FWFlux' : 0.0, 'HFlux':'26-50','Error':0.15} ,
             'Barents1' : {'Inflow' : -1.2, 'Outflow' : 4.1, 'Total flow' :'2-2.3', 'FWFlux' : 0.0, 'HFlux':'26-50','Error':0.15} ,
             'Barents2' : {'Inflow' : -1.2, 'Outflow' : 4.1, 'Total flow' :'2-2.3', 'FWFlux' : 0.0, 'HFlux':'26-50','Error':0.15} ,
             'Bering': {'Inflow' : 0, 'Outflow' : 0, 'Total flow' :0.8, 'FWFlux' : 0.0, 'HFlux': 0.0,'Error':0.2} ,
             'Davis': {'Inflow' : 0, 'Outflow' : 0, 'Total flow' :-2.6, 'FWFlux' : 0.0, 'HFlux':19.0,'Error':1.0},
             'Davis2': {'Inflow' : 0, 'Outflow' : 0, 'Total flow' :-2.6, 'FWFlux' : 0.0, 'HFlux':19.0,'Error':1.0},
             'Davis1': {'Inflow' : 0, 'Outflow' : 0, 'Total flow' :-2.6, 'FWFlux' : 0.0, 'HFlux':19.0,'Error':1.0}
            }

    obs.title = "Observ"
    obs.color = "grey"
    return liq10,core2,obs
