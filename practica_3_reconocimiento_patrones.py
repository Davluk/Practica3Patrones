#####################################################
#   bibliotecas necesarias para correr el codigo    #
#####################################################
import matplotlib.pyplot as plt
import random
import math

#####################################################
#        inialize the points to use                 #
#####################################################
print("init_program")
sample_points = []
for i in range(0,100000):
    sample_points.append([0,0])
    sample_points[i][0]=random.uniform(0,100)
    sample_points[i][1]=random.uniform(0,100)

#####################################################
#        functions to unse in the algorithm         #
#####################################################
def get_average_point(points):
    x_center = 0
    y_center = 0
    for unique_point in points:
        x_center+= unique_point[0]
        y_center+= unique_point[1]
    x_center/=(len(points)+0.0000001)
    y_center/=(len(points)+0.0000001)
    return [x_center,y_center]

def euc_dist(point_a,point_b):
    return math.sqrt( math.pow( point_b[0]-point_a[0] , 2 ) + math.pow( point_b[1]-point_a[1] , 2 ) )

def get_min_center_index(point,centers):
    min_index= 0
    min_dist = 100000000
    for index in range(0,len(centers)):
        l_dist = euc_dist( point , centers[index] )
        if( l_dist < min_dist ):
            min_index = index
            min_dist = l_dist
    return min_index
            
def dist_prom_area(center,area):
    distance = 0
    for point in area:
        distance+= euc_dist( point , center )
    #return distance/(len(area)+0.0000001)
    return distance

######################################################
#           LGB algorithm                            #
######################################################
def LGB(points,num_centers,E1,E2,tolerance):
    centers=[]
    areas=[]
    error = 0
    centers.append( get_average_point(points) )
    scaller = "\n"
    for index in range(100):
        scaller=scaller+"-"
    print(scaller)
    while len(centers)<num_centers:
        temp_centers = []
        gdist_ant = 0
        error = 10000000
        for center in centers:
            temp_centers.append( [ center[0]+E1[0] , center[1]+E1[1] ] )
            temp_centers.append( [ center[0]+E2[0] , center[1]+E2[1] ] )
        print()
        while error>tolerance:
            print(".")
            l_areas = []
            new_centers=[]
            temp_dist_average = 0
            #funciona
            for index in range(0,len(temp_centers)):
                l_areas.append([])

            for index in range(0,len(points)):
                l_areas[ get_min_center_index( points[index] , temp_centers ) ].append( points[index] )

            #funciona
            for index in range(0,len(l_areas)):
                new_centers.append( get_average_point( l_areas[index] ) )

            l_areas = []
            for index in range(0,len(temp_centers)):
                l_areas.append([])
            
            for index in range(0,len(points)):
                l_areas[ get_min_center_index( points[index] , new_centers ) ].append( points[index] )
            
            #funciona
            for index in range(0,len(new_centers)):
                temp_dist_average+=dist_prom_area(new_centers[index],l_areas[index])
            
            #error = math.fabs( temp_dist_average/len(new_centers) - gdist_ant )/(gdist_ant+1)
            error = math.fabs( temp_dist_average - gdist_ant )/(gdist_ant+1)
            #gdist_ant = temp_dist_average/len(new_centers)
            gdist_ant = temp_dist_average
            temp_centers = new_centers
            areas = l_areas

        centers = temp_centers
    return areas , centers

###########################################################
#           K-means algorithm                             #
###########################################################
def kmeans(points,num_centers,tolerance):
    centers=[]
    areas=[]
    error = 0
    for index in range(num_centers):
        centers.append( [random.random(),random.random()] )
    scaller = ""
    for index in range(100):
        scaller=scaller+"-"
    print(scaller)

    temp_centers = []
    gdist_ant = 0
    error = 100
    for center in centers:
        temp_centers.append( [ center[0] , center[1] ] )
        temp_centers.append( [ center[0] , center[1] ] )

    print()
    while error>tolerance:
        print(".",end='')
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
            new_centers.append( get_average_point( l_areas[index] ) )

        l_areas = []
        for index in range(0,len(temp_centers)):
            l_areas.append([])
        
        for index in range(0,len(points)):
            l_areas[ get_min_center_index( points[index] , new_centers ) ].append( points[index] )
        
        #funciona
        for index in range(0,len(new_centers)):
            temp_dist_average+=dist_prom_area(new_centers[index],l_areas[index])
        
        error = math.fabs( temp_dist_average - gdist_ant )/(gdist_ant+1)
        gdist_ant = temp_dist_average
        temp_centers = new_centers
        areas = l_areas

    centers = temp_centers
    return areas , centers


###########################################################
#           practica 3 - punto 1 a                        #
###########################################################
#print("init_lgb")

#for size in [4,8,16,64,256]:
#    lgb_areas ,lgb_centers = LGB( sample_points , size , [0.0,0.001] , [0.0,-0.001] , 0.01 )
#    X1 = []
#    Y1 = []
#    plt_areas = []

#    for center in lgb_centers:
#        X1.append(center[0])
#        Y1.append(center[1])

#    for index in range(len(lgb_areas)):
#        plt_areas.append([[],[]])
#        for point in lgb_areas[index]:
#            plt_areas[index][0].append( point[0] )
#            plt_areas[index][1].append( point[1] )

#    plt.figure()
#    ax1 = plt.subplot(111)
#    ax1.set_title('LGB '+ str(size))

#    for index in range(len(lgb_areas)):
#        ax1.scatter( plt_areas[index][0] , plt_areas[index][1] , alpha = 0.2)

#    ax1.scatter( X1 , Y1  ,c = 'black',alpha = 0.9)

#    ax1.axis([0 , 100 , 0 , 100])
#    ax1.grid(True)
#    plt.show()

###########################################################
#        practica 3 punto 1 b                             #
###########################################################
#print("init_lgb")
#for size in [4,8,16,64,256]:
#    lgb_areas ,lgb_centers = LGB( sample_points , size ,  [random.random(),random.random()] , [random.random(),random.random()] , 0.01 )
#    X1 = []
#    Y1 = []
#    plt_areas = []

#    for center in lgb_centers:
#        X1.append(center[0])
#        Y1.append(center[1])

#    for index in range(len(lgb_areas)):
#        plt_areas.append([[],[]])
#        for point in lgb_areas[index]:
#            plt_areas[index][0].append( point[0] )
#            plt_areas[index][1].append( point[1] )

#    plt.figure()
#    ax1 = plt.subplot(111)
#    ax1.set_title('LGB '+str(size))

#    for index in range(len(lgb_areas)):
#        ax1.scatter( plt_areas[index][0] , plt_areas[index][1] , alpha = 0.5)

#    ax1.scatter( X1 , Y1  ,c = 'white',alpha = 0.9)

#    ax1.axis([0 , 100 , 0 , 100])
#    ax1.grid(True)
#    plt.show()

###########################################################
#       practica 3 punto 1 - c                            #
###########################################################
print("init_lgb")
lgb_areas ,lgb_centers = kmeans( sample_points , 16 , 0.01 )
X1 = []
Y1 = []
plt_areas = []

for center in lgb_centers:
    X1.append(center[0])
    Y1.append(center[1])

for index in range(len(lgb_areas)):
    plt_areas.append([[],[]])
    for point in lgb_areas[index]:
        plt_areas[index][0].append( point[0] )
        plt_areas[index][1].append( point[1] )

plt.figure()
ax1 = plt.subplot(111)
ax1.set_title('Kmeans')

for index in range(len(lgb_areas)):
    ax1.scatter( plt_areas[index][0] , plt_areas[index][1] , alpha = 0.5)

ax1.scatter( X1 , Y1  ,c = 'white',alpha = 0.9)

ax1.axis([0 , 100 , 0 , 100])
ax1.grid(True)
plt.show()