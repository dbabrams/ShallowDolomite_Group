The Excel File contains pumpage data with the following fields:
1) p_num: unique identifier of each well
2) isws_facility_id: unique identifier of each facility
3) owner: name of the facility
4) fac_well_num: local number of the well at the facility
5) depth_total_last_known: the depth in feet of the well. Pay attention to this in case some wells have snuck in that are deeper than the aquifer you are modeling (hint hint, they did).
6 & 7) lam_x, lam_y: the x and y lambert coordinates. More on this after spring break.
8) remaining columns: the year with associated pumpage (in gallons reported for the year)

Recall that you need to fill in data gaps for years that a FACILITY (not just a well) did not report.
Recall that you need to remove OBVIOUS outliers.
You will need to spatially plot data for two years. Don't do this until you know which years you will be calibrating the model to (after spring break). 
