# This module will QAQC the current real-time data value.

"""
 Tests derived from QARTOD  Manual for Real-Time Quality Control of In-situ Temperature and Salinity Data
 Citation: U.S. Integrated Ocean Observing System, 2013. Manual for Real-Time Quality Control of In-situ Temperature and Salinity Data:
     A Guide to Quality Control and Quality Assurance of In-situ Temperature and Salinity Observations. 53pp.
"""
# Created by: Jason Adelaars (Moss Landing Marine Labs)
# Version: September 16, 2014

# Created for CeNCOOS weather and seawater data feeds with the goal of being able to QAQC mulitple parameters with this module


from datetime import datetime as dt
import numpy as np
import math


"""
The gap test determines that the most recent data point has been measured and received within
the expected time window (tim_inc) and has the correct time stamp (tim_stmp). This can only be used on real-time data
Cannot be used on post processing of historical data
Time format in UTC ISO 8601 (YYYY-MM-DDThh:mm:ssZ)
"""
def gapTest( tim_stmp , tim_inc ):
    # tim_stmp = datapoint's time stamp in UTC ISO 8601 format
    # tim_inc = sample frequency (in minutes)
    
    now = dt.utcnow()
    try:
        tim_stmp = dt.strptime(tim_stmp, '%Y-%m-%dT%H:%M:%SZ')
    except ValueError:
        raise ValueError("Timestamp '"+tim_stmp+"' not in correct format. Please use ISO 8601 (YYYY-MM-DDThh:mm:ssZ)")
    timedelta = now - tim_stmp
    mins = timedelta.seconds/60
    if mins > tim_inc:
        flag = 4
    else:
        flag = 1
    return flag

"""
Received data message contains the proper structure without any indicators of flawed
transmission such as parity errors. Possible tests are: a) the expected number of characters
for fixed length messages equals the number of characters received.
Best for CTD hex data strings
"""
def syntax( rec_char , nchar ):
    # rec_char = received datastring
    # nchar = expected number of characters
    if len(rec_char) != nchar:
        flag = 4
    else:
        flag = 1
    return flag

"""
Sensor range test evaluates the scaled sensor output value against the sensor's
minimum and maximum range
"""
def sensorRange( val , t_sensor_min , t_sensor_max ):
    # val = current data point
    # t_sensor_min = minimum scaled sensor value
    # t_sensor_max = maximum scaled sensor value
    if val < t_sensor_min or val > t_sensor_max:
        flag = 4
    elif math.isnan(val) == True:
        flag = 9
    else:
        flag = 1
    return flag


"""
User range test evaluates the scaled sensor output value against the a min/max range
defined by the user
"""
def userRange( val , t_user_min , t_user_max ):
    # val = current data point
    # t_user_min = minimum scaled sensor value
    # t_user_max = maximum scaled sensor value
    if val < t_user_min or val > t_user_max:
        flag = 3
    elif math.isnan(val) == True:
        flag = 9
    else:
        flag = 1
    return flag

"""
This check is for single value spikes, specifically the value at point n-1. Spikes consisting
of more than one datapoint are difficult to capture, but their onset may be flagged by the rate of change test.
The spike test consists of two operator-selected thresholds. Adjacent data points are averaged to form a spike reference.
The absolute value of the spike is tested to capture positive and negative spikes. Large spikes are easier to
identify and flag as failures. Smaller spikes may be real and are only flagged suspect. The thresholds may be fixed values or dynamic.

"""

def spike_ref( compList , thrshld_low , thrshld_high ):
    # compList = list of current plus previous two samples
    # thrshld_low = min user defined threshold
    # thrshld_high = max user defined threshold
    if len(compList) > 3:
        compList=compList[-3:]
    elif len(compList) < 3:
        flag = 9
        return flag
    mean = np.nanmean([compList[0],compList[2]])
    spk_ref = abs(compList[1]-mean)
    if spk_ref >= thrshld_high:
        flag = 4
    elif spk_ref < thrshld_high and spk_ref >= thrshld_low:
        flag = 3
    else:
        flag = 1

    return flag


##def spike_std( comp , thrshld_low=0 , thrshld_high=0 , dev=0 ):
##    # val = current data point
##    # comp = list of values to compare
##    # thrshld_low = min user defined threshold
##    # thrshld_high = max user defined threshold
##    # dev = number of deviations away from the mean user wants to define
##    # threshold defaults are 0 and test therefore defaults to dynamic thresholds if user thresholds are not set
##
##    if (thrshld_low == 0) and (thrshld_high == 0):
##        # if defaults used
##        if dev == 0:
##            raise ValueError("Please input a deviation multiplier (dev) for spike test.")
##        median = np.nanmedian(comp)
##        std = np.nanstd(comp)
##
##        if val > (median+(dev*std)) or val < (median-(dev*std)):
##            flag = 3
##        else:
##            flag = 1
##
##    else:
##        #if user-defined thresholds defined
##
##        mean = np.nanmean([comp[-2],val])
##        spk_ref = abs(comp[-1]-mean)
##        if spk_ref >= thrshld_high:
##            flag = 4
##        elif spk_ref < thrshld_high and spk_ref >= thrshld_low:
##            flag = 3
##        else:
##            flag = 1
##
##    return flag

"""
When some sensors and/or data collection platforms fail, the result can be a continuously repeated observation of the same value.
This test compares the present data point (n) to a number of previous observations. The current
datapoint is flagged if it has the same value as previous observations within a tolerance value, eps,
to allow for numerical round-off error.
"""
def flatLine( val , comp , eps , rep_cnt_fail=6 , rep_cnt_suspect=3 ):
    # val = current datapoint
    # comp = comparison values (n). Must be n>rep_cnt_fail
    # eps = tolerance value
    # rep_cnt_fail = number of repeated equal values to fail (default=6)
    # rep_cnt_suspect = number of repeated equal values to suspect (default=3)
    fail=rep_cnt_fail
    suspect=rep_cnt_suspect
    eps=eps
    if len(comp) < fail:
        raise ValueError("Not enough values to compare. Number of values in comp must be greater than rep_cnt_fail and rep_cnt_suspect.")
    

    count = 0
    flag=[]
    for i in np.arange(len(comp)):
        test = abs(val - comp[i])
        if test<=eps:
            count+= 1
    if count >= suspect:
        flag.append(3)

    count = 0
    for i in np.arange(len(comp)):
        test = abs(val - comp[i])
        if test<=eps:
            count+= 1
    if count >= fail:
        flag.append(4)

    else:
        flag.append(1)

    return max(flag)
"""
A common sensor failure mode can provide a data series that is nearly but not exactly a flat line
(e.g., if the sensor head were to become wrapped in debris). This test inspects for an SD value or a range variation
(MAX-MIN) value that fails to exceed threshold values (MIN_VAR_WARN, MIN_VAR_FAIL) over a user defined time period.
"""
def attenSig( comp , min_var_warn=0 , min_var_fail=0 , sd = False ):
    # comp: list of values to run the test on. Therefore, time interval determined by user
    # min_var_warn: minimum range of data variation to trip a warning
    # min_var_fail: minimum range of data variation to trip a fail
    # sd: if true, comp data stdev will be tested against stdev thresholds set by user
    if min_var_warn==0 and min_var_fail==0:
        raise ValueError ("Please indicate min&max deviance threshold values.")
    n=len(comp)

    if sd==True:
        if np.nanstd(comp) <= min_var_warn:
            flag = 3
        elif np.nanstd(comp) <= min_var_fail:
            flag = 4
        else:
            flag = 1


    if sd==False:
        test = abs(max(comp) - min(comp))
        if test <= min_var_warn:
            flag = 3
        elif test <= min_var_fail:
            flag = 4
        else:
            flag = 1

    flag_list = [flag*(x/x) for x in np.arange(1,(len(comp)+1))]

    return flag_list


    
        



    
    

        
            
        
        
    
    











