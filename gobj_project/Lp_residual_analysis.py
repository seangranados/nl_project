#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from analyze_residuals import Lp_analyze
import time 

# specify lp epochs to work on
epoch_lp_list = ["20040726nirc2", "20050730nirc2", "20050731nirc2",
             "20060502nirc2", "20090722nirc2", "20120516nirc2", "20130424nirc2", "20130425nirc2",
             "20140319nirc2", "20140320nirc2", "20140418nirc2", "20140511nirc2", "20140703nirc2",
             "20140803nirc2", "20140804nirc2", "20150331nirc2", "20160517nirc2", "20170716nirc2",
             "20170726nirc2", "20190514nirc2", "20190619nirc2", "20190814nirc2", "20200731nirc2",
             "20210713nirc2", "20210815nirc2", "20220515nirc2", "20230727nirc2"]

epoch_kp_list = ["20050731nirc2","20050731nirc2","20050731nirc2",
                 "20060502nirc2", "20090722nirc2", "20120515nirc2", "20130426_0427nirc2", "20130426_0427nirc2", 
                 "20140319nirc2", "20140320nirc2", "20140418nirc2", "20140512nirc2", "20140703nirc2", 
                 "20140803nirc2", "20140804nirc2","20150331nirc2", "20160503nirc2", "20170718nirc2", 
                 "20170727nirc2", "20190513nirc2", "20190625nirc2", "20190814nirc2", "20200809nirc2", 
                 "20210713nirc2", "20210815nirc2", "20220515nirc2","20230727nirc2"]


# specify matching or close kp epochs to work on

# got exception for align run for 0,1,2 align error

# 5,7 set of lp and kp

# create lp_res object, corr here is the correlation threshold for starfinder force mode when looking 
# for missed kp sources and kp/lp sources

#ran 16 already
    
# Manually process each index from 4 to 26, skipping 16

# lp_res5 = Lp_analyze(epoch_lp=epoch_lp_list[5], epoch_kp=epoch_kp_list[5], corr="0.4")
# lp_res5.sgra()
# lp_res5.setup_align()
# lp_res5.run_align()
# lp_res5.setup_force()
# lp_res5.missed_kp()
# lp_res5.align_force_compare()
# lp_res5.kp_and_lp_sources()
# lp_res5.assemble_final_list()

# lp_res6 = Lp_analyze(epoch_lp=epoch_lp_list[6], epoch_kp=epoch_kp_list[6], corr="0.4")
# lp_res6.sgra()
# lp_res6.setup_align()
# lp_res6.run_align()
# lp_res6.setup_force()
# lp_res6.missed_kp()
# lp_res6.align_force_compare()
# lp_res6.kp_and_lp_sources()
# lp_res6.assemble_final_list()

# lp_res7 = Lp_analyze(epoch_lp=epoch_lp_list[7], epoch_kp=epoch_kp_list[7], corr="0.4")
# lp_res7.sgra()
# lp_res7.setup_align()
# lp_res7.run_align()
# lp_res7.setup_force()
# lp_res7.missed_kp()
# lp_res7.align_force_compare()
# lp_res7.kp_and_lp_sources()
# lp_res7.assemble_final_list()

# lp_res8 = Lp_analyze(epoch_lp=epoch_lp_list[8], epoch_kp=epoch_kp_list[8], corr="0.4")
# lp_res8.sgra()
# lp_res8.setup_align()
# lp_res8.run_align()
# lp_res8.setup_force()
# lp_res8.missed_kp()
# lp_res8.align_force_compare()
# lp_res8.kp_and_lp_sources()
# lp_res8.assemble_final_list()

# lp_res9 = Lp_analyze(epoch_lp=epoch_lp_list[9], epoch_kp=epoch_kp_list[9], corr="0.4")
# lp_res9.sgra()
# lp_res9.setup_align()
# lp_res9.run_align()
# lp_res9.setup_force()
# lp_res9.missed_kp()
# lp_res9.align_force_compare()
# lp_res9.kp_and_lp_sources()
# lp_res9.assemble_final_list()

# lp_res10 = Lp_analyze(epoch_lp=epoch_lp_list[10], epoch_kp=epoch_kp_list[10], corr="0.4")
# lp_res10.sgra()
# lp_res10.setup_align()
# lp_res10.run_align()
# lp_res10.setup_force()
# lp_res10.missed_kp()
# lp_res10.align_force_compare()
# lp_res10.kp_and_lp_sources()
# lp_res10.assemble_final_list()

# lp_res11 = Lp_analyze(epoch_lp=epoch_lp_list[11], epoch_kp=epoch_kp_list[11], corr="0.4")
# lp_res11.sgra()
# lp_res11.setup_align()
# lp_res11.run_align()
# lp_res11.setup_force()
# lp_res11.missed_kp()
# lp_res11.align_force_compare()
# lp_res11.kp_and_lp_sources()
# lp_res11.assemble_final_list()

# lp_res12 = Lp_analyze(epoch_lp=epoch_lp_list[12], epoch_kp=epoch_kp_list[12], corr="0.4")
# lp_res12.sgra()
# lp_res12.setup_align()
# lp_res12.run_align()
# lp_res12.setup_force()
# lp_res12.missed_kp()
# lp_res12.align_force_compare()
# lp_res12.kp_and_lp_sources()
# lp_res12.assemble_final_list()

# lp_res13 = Lp_analyze(epoch_lp=epoch_lp_list[13], epoch_kp=epoch_kp_list[13], corr="0.4")
# lp_res13.sgra()
# lp_res13.setup_align()
# lp_res13.run_align()
# lp_res13.setup_force()
# lp_res13.missed_kp()
# lp_res13.align_force_compare()
# lp_res13.kp_and_lp_sources()
# lp_res13.assemble_final_list()

# lp_res14 = Lp_analyze(epoch_lp=epoch_lp_list[14], epoch_kp=epoch_kp_list[14], corr="0.4")
# lp_res14.sgra()
# lp_res14.setup_align()
# lp_res14.run_align()
# lp_res14.setup_force()
# lp_res14.missed_kp()
# lp_res14.align_force_compare()
# lp_res14.kp_and_lp_sources()
# lp_res14.assemble_final_list()

# lp_res15 = Lp_analyze(epoch_lp=epoch_lp_list[15], epoch_kp=epoch_kp_list[15], corr="0.4")
# lp_res15.sgra()
# lp_res15.setup_align()
# lp_res15.run_align()
# lp_res15.setup_force()
# lp_res15.missed_kp()
# lp_res15.align_force_compare()
# lp_res15.kp_and_lp_sources()
# lp_res15.assemble_final_list()

# # Skip 16

# lp_res17 = Lp_analyze(epoch_lp=epoch_lp_list[17], epoch_kp=epoch_kp_list[17], corr="0.4")
# lp_res17.sgra()
# lp_res17.setup_align()
# lp_res17.run_align()
# lp_res17.setup_force()
# lp_res17.missed_kp()
# lp_res17.align_force_compare()
# lp_res17.kp_and_lp_sources()
# lp_res17.assemble_final_list()

# lp_res18 = Lp_analyze(epoch_lp=epoch_lp_list[18], epoch_kp=epoch_kp_list[18], corr="0.4")
# lp_res18.sgra()
# lp_res18.setup_align()
# lp_res18.run_align()
# lp_res18.setup_force()
# lp_res18.missed_kp()
# lp_res18.align_force_compare()
# lp_res18.kp_and_lp_sources()
# lp_res18.assemble_final_list()

# lp_res19 = Lp_analyze(epoch_lp=epoch_lp_list[19], epoch_kp=epoch_kp_list[19], corr="0.4")
# lp_res19.sgra()
# lp_res19.setup_align()
# lp_res19.run_align()
# lp_res19.setup_force()
# lp_res19.missed_kp()
# lp_res19.align_force_compare()
# lp_res19.kp_and_lp_sources()
# lp_res19.assemble_final_list()

# lp_res20 = Lp_analyze(epoch_lp=epoch_lp_list[20], epoch_kp=epoch_kp_list[20], corr="0.4")
# lp_res20.sgra()
# lp_res20.setup_align()
# lp_res20.run_align()
# lp_res20.setup_force()
# lp_res20.missed_kp()
# lp_res20.align_force_compare()
# lp_res20.kp_and_lp_sources()
# lp_res20.assemble_final_list()

# lp_res21 = Lp_analyze(epoch_lp=epoch_lp_list[21], epoch_kp=epoch_kp_list[21], corr="0.4")
# lp_res21.sgra()
# lp_res21.setup_align()
# lp_res21.run_align()
# lp_res21.setup_force()
# lp_res21.missed_kp()
# lp_res21.align_force_compare()
# lp_res21.kp_and_lp_sources()
# lp_res21.assemble_final_list()

# lp_res22 = Lp_analyze(epoch_lp=epoch_lp_list[22], epoch_kp=epoch_kp_list[22], corr="0.4")
# lp_res22.sgra()
# lp_res22.setup_align()
# lp_res22.run_align()
# lp_res22.setup_force()
# lp_res22.missed_kp()
# lp_res22.align_force_compare()
# lp_res22.kp_and_lp_sources()
# lp_res22.assemble_final_list()

# lp_res23 = Lp_analyze(epoch_lp=epoch_lp_list[23], epoch_kp=epoch_kp_list[23], corr="0.4")
# lp_res23.sgra()
# lp_res23.setup_align()
# lp_res23.run_align()
# lp_res23.setup_force()
# lp_res23.missed_kp()
# lp_res23.align_force_compare()
# lp_res23.kp_and_lp_sources()
# lp_res23.assemble_final_list()

# lp_res24 = Lp_analyze(epoch_lp=epoch_lp_list[24], epoch_kp=epoch_kp_list[24], corr="0.4")
# lp_res24.sgra()
# lp_res24.setup_align()
# lp_res24.run_align()
# lp_res24.setup_force()
# lp_res24.missed_kp()
# lp_res24.align_force_compare()
# lp_res24.kp_and_lp_sources()
# lp_res24.assemble_final_list()

# lp_res25 = Lp_analyze(epoch_lp=epoch_lp_list[25], epoch_kp=epoch_kp_list[25], corr="0.4")
# lp_res25.sgra()
# lp_res25.setup_align()
# lp_res25.run_align()
# lp_res25.setup_force()
# lp_res25.missed_kp()
# lp_res25.align_force_compare()
# lp_res25.kp_and_lp_sources()
# lp_res25.assemble_final_list()

lp_res26 = Lp_analyze(epoch_lp=epoch_lp_list[26], epoch_kp=epoch_kp_list[26], corr="0.4")
lp_res26.sgra()
lp_res26.setup_align()
lp_res26.run_align()
lp_res26.setup_force()
lp_res26.missed_kp()
lp_res26.align_force_compare()
lp_res26.kp_and_lp_sources()
lp_res26.assemble_final_list()




# lp_res.starplant()
# lp_res.create_residual()
