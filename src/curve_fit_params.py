#curve fitting helper functions with parameters for PacBio and IonTorrent Datasets 
import numpy 
import random
import scipy
import sys
from numpy import array
from scipy import interpolate
from scipy.optimize import *
import math
numpy.seterr(over='ignore')
''' basic function forms for curve fitting'''

def linear(x,m,b):
    if len(x)==1: 
        return numpy.resize(m*numpy.array(x)+b,2)
    else:
        return m*numpy.array(x)+b 

def poly2(vals,a,b,c):
    output = [] 
    for x in vals:
        output.append(a*math.pow(x,2)+b*x+c)
    return output

def poly3(vals,a,b,c,d):
    output = [] 
    for x in vals: 
        output.append(a*math.pow(x,3)+b*pow(x,2)+c*x+d)
    return output 

def poly4(vals,a,b,c,d,e):
    output = [] 
    for x in vals:
        output.append(a*math.pow(x,4)+b*pow(x,3)+c*pow(x,2)+d*x+e)
    return output 

def poly5(vals,a,b,c,d,e,f):
    output =[] 
    for x in vals: 
        output.append(a*math.pow(x,5)+b*pow(x,4)+c*pow(x,3)+d*pow(x,2)+e*x+f)
    return output 

def exp(vals,a,b):
    output=[] 
    for x in vals: 
        output.append(a*numpy.exp(b*numpy.array(x)))
    return output 

def exp_sum(vals,a,b,c,d):
    output = [] 
    for x in vals: 
        output.append(a*numpy.exp(b*numpy.array(x))+c*numpy.exp(d*numpy.array(x)))
    return output

def log(vals,a,b):
    output=[] 
    for x in vals: 
        if (numpy.array(x)*b) < 0.001: 
            output.append(a*numpy.log(0.001)); 
        else:
            output.append(a*numpy.log(numpy.array(x)*b)) 
    return output 

def power(vals,a,b):
    output = [] 
    for x in vals: 
        output.append(a*numpy.power(numpy.array(x),b))
    return output 

def gauss(x, A, mu,sigma):
    return numpy.array(A*numpy.exp(-(x-mu)**2/(2.*sigma**2)))



#since it is unclear which type of function gives best fit, try common types and choose the one that gives 
#the least error. 
def fitToData(data,rangeToUse,gauss_params,numblocks,datatype):
    #split the data up into a range for empirical fit and for extrapolation. 
    all_keys = [] 
    all_fitted = []
    keys_used=[]
    data_fitted=[]
    datamax = int(max(data.keys()))
    for block in range(numblocks):
        #range for empirical data 
        start = int(rangeToUse[1]/float(numblocks))*block
        cap = int(int(rangeToUse[1])/float(numblocks)*(block+1))
        emp_keys=[] 
        emp_vals =[]

        for entry in data: 
            if entry > start and entry < cap+1: 
                emp_keys.append(entry)
                emp_vals.append(data[entry])
        emp_keys= numpy.array(emp_keys) 
        emp_vals = numpy.array(emp_vals) 
        #fit  common function forms to the data, determine which function form minimizes the least squares error. 
        try:
            popt_linear, pcov_linear = curve_fit(linear, emp_keys, emp_vals)
            lin_vals=numpy.array(linear(emp_keys,popt_linear[0],popt_linear[1]))
            lin_error = sum(abs(lin_vals-emp_vals)) 
            min_err = lin_error 
        #interpolate missing values and extrapolate to remaining values 
            if block == numblocks-1: 
                data_fitted = linear([i for i in range(start,datamax)],popt_linear[0],popt_linear[1])
                keys_used=[i for i in range(start,datamax)]
                
            else:
                data_fitted = linear([i for i in range(start,cap+1)],popt_linear[0],popt_linear[1])
                keys_used= [i for i in range(start,cap+1)]
        except Exception as e:
            #print " linear fit gave error:"+str(e)+'\n' 
            failure='linear fit did not work.\n'
        try:
            popt_poly2,pcov_poly2 = curve_fit(poly2,emp_keys,emp_vals)
            poly2_vals = numpy.array(poly2(emp_keys,popt_poly2[0],popt_poly2[1],popt_poly2[2]))
            poly2_error = sum(abs(poly2_vals - emp_vals)) 
            if poly2_error < min_err: 
                min_err = poly2_error
                if block == numblocks-1:
                    data_fitted=poly2([i for i in range(start,datamax)],popt_poly2[0],popt_poly2[1],popt_poly2[2])
                    keys_used=[i for i in range(start,datamax)]
                else: 
                    data_fitted=poly2([i for i in range(start,cap+1)],popt_poly2[0],popt_poly2[1],popt_poly2[2])
                    keys_used=[i for i in range(start,cap+1)]
        except Exception as e:
            #print "poly2 fit gave error:"+str(e)+'\n' 
            failure= "poly2 fit didn't work\n"
        try:
            popt_poly3,pcov_poly3 = curve_fit(poly3,emp_keys,emp_vals)
            poly3_vals=numpy.array(poly3(emp_keys,popt_poly3[0],popt_poly3[1],popt_poly3[2],popt_poly3[3]))
            poly3_error = sum(abs(poly3_vals - emp_vals)) 
            #print "poly3 error: " + str(poly3_error)+"\n"
            if poly3_error < min_err: 
                min_err = poly3_error 
                if block ==numblocks-1:
                    data_fitted=poly3([i for i in range(start,datamax)],popt_poly3[0],popt_poly3[1],popt_poly3[2],popt_poly3[3])
                    keys_used=[i for i in range(start,datamax)]
                else:
                    data_fitted=poly3([i for i in range(start,cap+1)],popt_poly3[0],popt_poly3[1],popt_poly3[2],popt_poly3[3])
                    keys_used=[i for i in range(start,cap+1)]
        except Exception as e:
            #print "poly3 fit gave error:"+str(e)+'\n' 
            failure="poly3 fit didn't work\n"
        try:
            popt_poly4,pcov_poly4 = curve_fit(poly4,emp_keys,emp_vals)
            poly4_vals=numpy.array(poly4(emp_keys,popt_poly4[0],popt_poly4[1],popt_poly4[2],popt_poly4[3],popt_poly4[4]))
            poly4_error = sum(abs(poly4_vals - emp_vals)) 
            if poly4_error < min_err: 
                min_err=poly4_error 
                if block == numblocks-1:              
                    data_fitted=poly4([i for i in range(start,datamax)],popt_poly4[0],popt_poly4[1],popt_poly4[2],popt_poly4[3],popt_poly4[4])
                    keys_used=[i for i in range(start,datamax)]
                else: 
                    data_fitted=poly4([i for i in range(start,cap+1)],popt_poly4[0],popt_poly4[1],popt_poly4[2],popt_poly4[3],popt_poly4[4])
                    keys_used= [i for i in range(start,cap+1)]
        except Exception as e:
            #print "poly4 fit gave error:"+str(e)+'\n'
            failure= "poly4 fit didn't work\n"
        try:
            popt_poly5,pcov_poly5 = curve_fit(poly5,emp_keys,emp_vals)
            poly5_vals=numpy.array(poly5(emp_keys,popt_poly5[0],popt_poly5[1],popt_poly5[2],popt_poly5[3],popt_poly5[4],popt_poly5[5]))
            poly5_error=sum(abs(poly5_vals - emp_vals)) 
            if poly5_error < min_err: 
                min_err = poly5_error 
                if block==numblocks-1: 
                    data_fitted=poly5([i for i in range(start,datamax)],popt_poly5[0],popt_poly5[1],popt_poly5[2],popt_poly5[3],popt_poly5[4],popt_poly5[5])
                    keys_used=[i for i in range(start,datamax)]
                else:   
                    data_fitted=poly5([i for i in range(start,cap+1)],popt_poly5[0],popt_poly5[1],popt_poly5[2],popt_poly5[3],popt_poly5[4],popt_poly5[5])
                    keys_used= [i for i in range(start,cap+1)]
        except Exception as e:
            #print "poly5 fit gave erorr:"+str(e)+'\n' 
            failure= "poly5 fit didnt' work\n"
            
        try:
            popt_exp,pcov_exp = curve_fit(exp,emp_keys,emp_vals)
            exp_vals = exp(emp_keys,popt_exp[0],popt_exp[1]) 
            exp_error=sum(abs(exp_vals-emp_vals)) 
            if exp_error < min_err: 
                min_err=exp_error
                if block == numblocks-1:
                    data_fitted= exp([i for i in range(start,datamax)],popt_exp[0],popt_exp[1]) 
                    keys_used=[i for i in range(start,datamax)]
                else:
                    data_fitted= exp([i for i in range(start,cap+1)],popt_exp[0],popt_exp[1])
                    keys_used=[i for i in range(start,cap+1)]
        except Exception as e:
            #print "exponential fit gave error:"+str(e)+'\n' 
            failure= "exponential fit didn't work\n"
        try:
            popt_exp_sum,pcov_exp_sum = curve_fit(exp_sum,emp_keys,emp_vals)
            exp_sum_vals = exp_sum(emp_keys,popt_exp_sum[0],popt_exp_sum[1],popt_exp_sum[2],popt_exp_sum[3]) 
            exp_sum_error=sum(abs(exp_sum_vals-emp_vals)) 
            if exp_sum_error < min_err: 
                min_err=exp_sum_error
                if block==numblocks-1:
                    data_fitted= exp_sum([i for i in range(start,datamax)],popt_exp_sum[0],popt_exp_sum[1],popt_exp_sum[2],popt_exp_sum[3]) 
                    keys_used=[i for i in range(start,datamax)]
                else:
                    data_fitted= exp_sum([i for i in range(start,cap+1)],popt_exp_sum[0],popt_exp_sum[1],popt_exp_sum[2],popt_exp_sum[3]) 
                    keys_used= [i for i in range(start,cap+1)]
        except Exception as e:
            #print "exponential sum fit gave error:"+str(e)+'\n' 
            failure= "exponential sum fit didn't work\n"
        
        try:
            popt_log,pcov_log = curve_fit(log,emp_keys,emp_vals)
            log_vals = log(emp_keys,popt_log[0],popt_log[1]) 
            log_error = sum(abs(log_vals-emp_vals)) 
            if log_error < min_err:
                min_err=log_error 
                if block==numblocks-1:    
                    data_fitted = log([i for i in range(start,datamax)],popt_log[0],popt_log[1]) 
                    keys_used=[i for i in range(start,datamax)]
                else:
                    data_fitted = log([i for i in range(start,cap+1)],popt_log[0],popt_log[1]) 
                    keys_used=[i for i in range(start,cap+1)]
        except Exception as e:
            #print "logarithmic fit gave error:"+str(e)+'\n' 
            failure= "logarithmic fit didn't work\n"
        try:
            popt_power,pcov_power=curve_fit(power,emp_keys,emp_vals)
            pow_vals = pow(emp_keys,popt_power[0],popt_power[1]) 
            pow_error= sum(abs(pow_vals-emp_vals)) 
            if pow_error < min_err: 
                min_err=pow_error
                if block==numblocks-1:
                    data_fitted = pow([i for i in range(start,datamax)],popt_power[0],popt_power[1]) 
                    keys_used= [i for i in range(start,datamax)]
                else:
                    data_fitted = pow([i for i in range(start,cap+1)],popt_power[0],popt_power[1]) 
                    keys_used= [i for i in range(start,cap+1)]
        except Exception as e:
            #print "power fit gave error:"+str(e)+'\n' 
            failure= "power fit didn't work\n"
        try: 
            popt_gauss,pcov_power = curve_fit(gauss,emp_keys,emp_vals,gauss_params)
            gauss_vals = gauss(emp_keys,popt_gauss[0],popt_gauss[1],popt_gauss[2])
            gauss_error = sum(abs(gauss_vals-emp_vals))
            if gauss_error < min_err: 
                min_err = gauss_error 
                if block==numblocks-1: 
                    data_fitted = gauss([i for i in range(start,datamax)],popt_gauss[0],popt_gauss[1],popt_gauss[2])
                    keys_used=[i for i in range(start,datamax)]
                else: 
                    data_fitted = gauss([i for i in range(start,cap+1)],popt_gauss[0],popt_gauss[1],popt_gauss[2])
                    keys_used =  [i for i in range(start,cap+1)]
        except Exception as e:
            #print "gaussian fit gave error:"+str(e)+'\n' 
            failure= "gaussian fit didn't work\n"
        
        all_keys = all_keys+list(keys_used)
        all_fitted=all_fitted+list(data_fitted) 
        #extrapolate from the final function 
    fitDict = dict() 
    last_non_zero_index=0 
    for i in range(len(all_keys)): 
        if all_fitted[i] < 0: 
            next_non_zero_index = None 
            for j in range(1,len(all_keys)-(i)):
                if all_fitted[i+j] > 0: 
                    next_non_zero_index = i+j
                    break 
            if next_non_zero_index ==None:
                next_non_zero_index = last_non_zero_index 
            fitDict[all_keys[i]]=.5*(all_fitted[last_non_zero_index]+all_fitted[next_non_zero_index])
        else:
            last_non_zero_index=i 
            fitDict[all_keys[i]] = all_fitted[i]
    return [fitDict,datatype]  
    
        
