#####################################################
#   bibliotecas necesarias para correr el codigo    #
#####################################################
from numpy.linalg import inv
import matplotlib.pyplot as plt
import sounddevice as sd
import scipy.io.wavfile as wav
import numpy as np
import random
from scipy import signal
import math
import csv
import sys
import time
print("init_program")

########################################################
#           constants of the program
########################################################
files = ['uno','dos','tres','cuatro','cinco','seis','siete','ocho','nueve','diez']
hamming_size = 160
hamming_jump = 64
r_vect_size = 12
umbral_1 = 0.000004
umbral_2 = 0.00000027

recs = []
preemph_sigs = [] 
signal_chunks= []
sig_power = []
signal_limits = []
true_signal=[]
signal_r_predict = []
a_vector_centers = []
quantizers_points = []
quantizers_centers = []
quantizers_areas = []

#########################################################
#   funcion para calcular el vector de correlaciones    #
#########################################################
def Rxx(x,i,M):                   #vector de muestra y tamaño de la ventana
    suma = 0
    for k in range(0,M-1):
        if((k+i)>0 and (k+i)<len(x)):
            suma += x[k]*x[k+i]/M
    return suma

def get_r_vector(x,maximum):
    temp_vector=[]
    for index in range(0,maximum):
        temp_vector.append( Rxx( x , index , len( x ) ) )
    return temp_vector

def get_R_vector(X,maximum):
    temp_vector = []
    for i in range( 0, maximum ):
        temp_vector.append([])
        for j in range( 0, maximum ):
            temp_vector[i].append( Rxx( X , int(math.fabs(i-j)) , len( X ) ) )
    return temp_vector

def get_W_vector(senal,maximum):
    senal_r_vector = []
    senal_R_vector = []
    senal_w_vector = []
    for index in range( 0 , len(senal) ):
        senal_r_vector.append( get_r_vector( senal[index] , maximum ) )
    for index in range( 0 , len(senal) ):
        senal_R_vector.append( get_R_vector( senal[index] , maximum ) )
    for index in range( 0 , len(senal) ):
        temp_r_array = np.array( senal_r_vector[index] )
        temp_R_array = np.array( senal_R_vector[index] )
        senal_w_vector.append( np.dot( inv( temp_R_array ) , np.transpose( temp_r_array ) ) )
    return senal_w_vector

def get_w_vector(senal,maximum):
    temp_r_vector = get_r_vector( senal , maximum )
    temp_R_vector = get_R_vector( senal , maximum )
    return np.dot( inv( temp_R_vector ) , np.transpose( temp_r_vector ) )

def get_w_from_r(r):
    R_matrix=[]
    for i_index in range(len(r)-1):
        R_matrix.append([])
        for j_index in range(len(r)-1):
            R_matrix[i_index].append( r[int(abs(j_index-i_index))] )
    temp_r = np.array( r[1:len(r)] )
    temp_R = np.array( R_matrix )
    return np.dot( inv( temp_R ) , np.transpose( temp_r ) ) 


def get_a_from_w(W):
    temp_a = []
    temp_a.append( 1 )
    for index in range(len(W)):
        temp_a.append( -W[index] )
    return temp_a


def get_init_end_signal(pwr_signal,umbral_1,umbral_2):
    init_signal=0
    end_signal=0
    for index in range(0,len(pwr_signal)):
        if(pwr_signal[index]>umbral_1 and init_signal==0 and index>50 and init_signal==0):
            init_signal = index
        #if((index+2)<len(pwr_signal) and pwr_signal[index]<umbral_2 and end_signal==0 and init_signal!=0 and pwr_signal[index+2]<umbral_2):
        if(init_signal!=0 and (index+2)<len(pwr_signal) and pwr_signal[index]>umbral_2 and pwr_signal[index+1]>umbral_2 and pwr_signal[index+2]>umbral_2 ):
            end_signal = index
    return init_signal,end_signal


def pre_enphasis_filter(X): 
    temp_vector=[]
    for i in range(0,len(X)):
        if((i-1)>0):
            temp_vector.append( X[i] - 0.95*X[i-1] )
        else:
            temp_vector.append( X[i] )
    return temp_vector

Hamming_window = np.hamming(hamming_size)
def get_vector_with_hamming(J,vector):
    temp_vector = []
    for i in range(0,hamming_size):
        if((J+i)<len(vector)):
            temp_vector.append( vector[J+i]*Hamming_window[i] )
        else:
            temp_vector.append( 0 )
    return temp_vector


def get_average_point(points):
    l_average_point=[]
    temp_center = 0
    for i_index in range(0,r_vect_size):
        temp_center = 0
        for p_index in range( len( points ) ):
            temp_center+= points[p_index][i_index]
        temp_center/=(len(points)+0.0000001)
        l_average_point.append(temp_center)
    return l_average_point

#it is doesnt used 
def euc_dist(point_a,point_b):
    suma = 0
    for index in range(r_vect_size):
        suma += math.pow(point_b[index]-point_a[index],2)
    return math.sqrt( suma )

def get_min_center_index(point,centers):
    min_index= 0
    min_dist = 100000000000000
    for index in range(0,len(centers)):
        if( get_itakura_saito_distance( centers[index] , point ) < min_dist ):
            min_index = index
            min_dist = get_itakura_saito_distance( centers[index] , point )
    return min_index

def dist_prom_area(center,area):
    distance = 0
    for point in area:
        distance+= get_itakura_saito_distance( center , point )
    return distance/(len(area)+0.000000000001)

def get_itakura_saito_distance( r_vector_model , r_vector_sample ):
    w_vector = get_w_from_r( r_vector_model )
    a_vector = get_a_from_w( w_vector )
    dist_int = Rxx( a_vector , 0 , len(a_vector) )*r_vector_sample[0]
    for index in range( 1 , r_vect_size ):
        dist_int+= 2*Rxx( a_vector , index , len(a_vector) )*r_vector_sample[index]
    return dist_int

def LGB(points,num_centers,E1,E2,tolerance):

    centers=[]
    areas=[]
    error = 0
    centers.append( get_average_point(points) )
    scaller = ""
    for index in range(100):
        scaller=scaller+"-"
    print(scaller)

    while len(centers)<num_centers:
        temp_centers = []
        gdist_ant = 0
        error = 100
        for center in centers:
            tmp_vector1=[]
            tmp_vector2=[]
            for i_index in range( r_vect_size ):
                tmp_vector1.append( center[i_index]+E1[i_index] )
                tmp_vector2.append( center[i_index]+E2[i_index] )
            temp_centers.append( tmp_vector1 )
            temp_centers.append( tmp_vector2 )
        #print(temp_centers," centers")
        while error>tolerance:
        #for count in range(30):
            l_areas = []
            new_centers=[]
            temp_dist_average = 0
            #funciona
            for index in range(0,len(temp_centers)):
                l_areas.append([])

            for index in range(0,len(points)):
                l_areas[ get_min_center_index( points[index] , temp_centers ) ].append( points[index] )

            #funciona

            for index in range(0,len(temp_centers)):
                if(len(l_areas[index])):
                    new_centers.append( get_average_point( l_areas[index] ) )

            l_areas = []
            for index in range(0,len(new_centers)):
                l_areas.append([])
            
            for index in range(0,len(points)):
                l_areas[ get_min_center_index( points[index] , new_centers ) ].append( points[index] )
            
            #funciona
            for index in range(0,len(new_centers)):
                temp_dist_average+=dist_prom_area(new_centers[index],l_areas[index])
            
            print()
            for index in range(len(l_areas)):
                print("- ",index,":",len(l_areas[index]),end=' ')

            error = math.fabs( temp_dist_average - gdist_ant )/(gdist_ant+1)
            gdist_ant = temp_dist_average
            temp_centers = new_centers
            areas = l_areas

        centers = temp_centers
    return areas , centers

##########################################################
#       practica 
##########################################################
t_start = time.time()
t_end = 0
for index_files in range(len(files)):
    recs.append([])
    for index in range(15):
        recs[index_files].append([])
        with open('recs/'+files[index_files]+'/'+files[index_files]+str(index)+'.csv','r') as f:
            reader = csv.reader(f)
            for row in reader:
                if(len(row)>0):
                    recs[index_files][index].append(float(row[0]))
    print(files[index_files]+"\t samples readed")
t_end = time.time()
print("execution time ",int((t_end-t_start)/3600)," hours",int((t_end-t_start)/60),' minutes',int((t_end-t_start)%60),' seconds' )

print("#############################################\nPRE-EMPHASIS")
t_start = time.time()
for index_recs in range(len(recs)):
    print(files[index_recs]+'\tsamples: applying pre-emphasis')
    preemph_sigs.append([])
    for index_samps in range(15):
        preemph_sigs[index_recs].append([])
        preemph_sigs[index_recs][index_samps] = pre_enphasis_filter( recs[index_recs][index_samps] )
t_end = time.time()
print("execution time ",int((t_end-t_start)/3600)," hours",int((t_end-t_start)/60),' minutes',int((t_end-t_start)%60),' seconds' )


print("\n############################################\nHAMMING-WINDOWS")
t_start = time.time()
jumps = int( len(recs[0][0] )/hamming_jump )
for rec_index in range(len(preemph_sigs)):
    print(files[rec_index]+'\tsamples: applying hamming windowing')
    signal_chunks.append([])
    for samp_index in range(len(preemph_sigs[rec_index])):
        signal_chunks[rec_index].append([])
        for jump in range(jumps):
            signal_chunks[rec_index][samp_index].append( get_vector_with_hamming( jump*hamming_jump , preemph_sigs[rec_index][samp_index] ) )
t_end = time.time()
print("execution time ",int((t_end-t_start)/3600)," hours",int((t_end-t_start)/60),' minutes',int((t_end-t_start)%60),' seconds' )

print("\n############################################\nPOWER OF THE SIGNAL")
t_start = time.time()
for rec_index in range( len( signal_chunks ) ):
    print(files[rec_index]+'\tsamples: getting que power of every window')
    sig_power.append([])
    #for samp_index in range( signal_chunks[rec_index] ):
    for samp_index in range( 10 ):
        sig_power[rec_index].append([])
        for jump in range( len(signal_chunks[rec_index][samp_index]) ):
            sig_power[rec_index][samp_index].append( Rxx( signal_chunks[rec_index][samp_index][jump] , 0 , r_vect_size ) )
t_end = time.time()
print("execution time ",int((t_end-t_start)/3600)," hours",int((t_end-t_start)/60),' minutes',int((t_end-t_start)%60),' seconds' )

print("\n############################################\nBEGIN END SIGNAL")
t_start = time.time()
for rec_index in range( len( sig_power ) ):
    print(files[rec_index]+'\tsamples: getting the limmits of every signal')
    signal_limits.append([])
    for samp_index in range( len(sig_power[rec_index]) ):
        signal_limits[rec_index].append([])
        temp_begin,temp_end = get_init_end_signal( sig_power[rec_index][samp_index] , umbral_1 , umbral_2 )
        #print("begin:\n",temp_begin,"\tend:\t",temp_end)
        signal_limits[rec_index][samp_index].append( temp_begin )
        signal_limits[rec_index][samp_index].append( temp_end )

t_end = time.time()
print("execution time ",int((t_end-t_start)/3600)," hours",int((t_end-t_start)/60),' minutes',int((t_end-t_start)%60),' seconds' )

print("\n############################################\nCUTTING")
t_start = time.time()

for rec_index in range( len( signal_chunks ) ):
    print(files[rec_index]+'\tsamples: getting the procesable signal ')
    true_signal.append([])
    for samp_index in range( 10 ):
        true_signal[rec_index].append([])
        true_signal[rec_index][samp_index]=signal_chunks[rec_index][samp_index][signal_limits[rec_index][samp_index][0]:signal_limits[rec_index][samp_index][1]]

t_end = time.time()
print("execution time ",int((t_end-t_start)/3600)," hours",int((t_end-t_start)/60),' minutes',int((t_end-t_start)%60),' seconds' )

print("\n############################################\nr VECTORS")
t_start = time.time()
for rec_index in range(len(true_signal)):
    print(files[rec_index]+'\tsamples: getting the r vectors for every chunk')
    signal_r_predict.append( [] )
    for samp_index in range(10):
        signal_r_predict[rec_index].append( [] )    
        for jump in range( len( true_signal[rec_index][samp_index] ) ):
            signal_r_predict[rec_index][samp_index].append( get_r_vector( true_signal[rec_index][samp_index][jump], r_vect_size+1 ) )
t_end = time.time()
print("execution time ",int((t_end-t_start)/3600)," hours",int((t_end-t_start)/60),' minutes',int((t_end-t_start)%60),' seconds' )


# concatente all the windows 
print("\n############################################\nSETTING THE WINDOWS in a bag")
signal_r_join = []
t_start = time.time()
for rec_index in range(len(true_signal)):
    print(files[rec_index]+'\tsamples: setting every window on a bag')
    signal_r_join.append( [] )
    for samp_index in range(10):
        for jump in range( signal_limits[rec_index][samp_index][1]-signal_limits[rec_index][samp_index][0] ):
            signal_r_join[rec_index].append( signal_r_predict[rec_index][samp_index][jump] )  
        
t_end = time.time()
print("execution time ",int((t_end-t_start)/3600)," hours",int((t_end-t_start)/60),' minutes',int((t_end-t_start)%60),' seconds' )

tmp_e1 = []
tmp_e2 = []
tmp_e1.append(  0.0001 ) 
tmp_e2.append( -0.0001 ) 
for i_index in range(1,r_vect_size):
        tmp_e1.append(  0.0001 ) 
        tmp_e2.append( -0.0001 )        
print("\n############################################\nCENTROIDS")
for rec_index in range( len( signal_r_join ) ):
    print(files[rec_index]+'\tsamples: getting the centroids and areas from every quantifier')
    temp_centers, temp_areas = LGB( signal_r_join[rec_index] , 8 , tmp_e1 , tmp_e2 , 0.000001 )
    quantizers_centers.append( temp_centers )
    quantizers_areas.append( temp_areas )





