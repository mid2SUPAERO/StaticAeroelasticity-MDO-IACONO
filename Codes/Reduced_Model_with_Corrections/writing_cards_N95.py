"""
NASTRAN CARDS WRITING
"""

#%reset to clear all variables
# clear to clear console

#print('\n'*50)
#print(" ======== Execution ======== ")
#print(" ")



# Possible issue: no more than 8 characters allowed in free field




# SOL card ====================================================================
def SOL(type):
    sol_card = ''
    if type=='static':
        sol_card += 'SOL,101\nCEND\n\n'
    elif type=='modal':
        sol_card += 'SOL,103\nCEND\n\n'
    else:
        print('wrong solution type\n\n')
        sol_card += 'fatal error'     
    return sol_card





# Boundary Condition: clamped end =============================================
def BC_clamped(node):
    BC1  = 'SPC = 1\n'
    BC2  = 'SPC1,1,123456,'+str(node)+'\n'
    #BC2 += 'SPCADD, 2, 1\n\n'
    return(BC1,BC2)



# MAT1 card ===================================================================
def MAT1(MID,E,G,rho):
    _ = ' ' #blank/empty
    mat1_card = 'MAT1'+','+str(MID)+','+str(E)+','+str(G)+','+_+','+str(rho)+'\n\n'
    return mat1_card





# This card has been modified in S3!
# GRID card ===================================================================
# SRef_Geo: SRef in which geometry is defined (grid point coordinates)
# SRef_Ana: SRef in which analyses are performed
def GRID(ID,SRef_Geo,X,SRef_Ana):
    _ = ' ' #blank/empty
    [X1,X2,X3]=[str(X[0]),str(X[1]),str(X[2])]
    grid_card = 'GRID'+','+str(ID)+','+ str(SRef_Geo) +','+X1+','+X2+','+X3+','+str(SRef_Ana)+'\n'
    return grid_card




# CBEAM card ==================================================================
# CBEAM, EID, PID,  GA,  GB,  X1,  X2,  X3,    ,
#      ,    ,    , W1A, W2A, W3A, W1B, W2B, W3B,
# EID      -> element id
# PID      -> property id
# GA,GB    -> the two grid points of the beam element
# X1,X2,X3 -> orientation vector of local y-axis (x=local beam axis)
# WA,WB    -> offset vectors from GA,GB of the shear center


def CBEAM(EID,PID,GA,GB,X,W):
    
    # inputs: 
    # EID,PID,GA,GB -> integers
    # X,WA,WB -> vectors (3 components)
             
    # Conversion to strings
    EID = str(EID)
    PID = str(PID)
    GA  = str(GA)
    GB  = str(GB)
    [X1,X2,X3]=[str(X[0]),str(X[1]),str(X[2])]
    [W1A,W2A,W3A]=[str(W[0]),str(W[1]),str(W[2])]
    [W1B,W2B,W3B]=[str(W[0]),str(W[1]),str(W[2])]
    _ = ' ' #blank/empty
       
    # Generation of the card
    cbeam_card  = ''
    cbeam_card += 'CBEAM'+','+EID+','+PID+','+GA +','+GB +','+X1 +','+X2 +','+X3 +','+ _ +', \n'\
               +  '     '+','+ _ +','+ _ +','+W1A+','+W2A+','+W3A+','+W1B+','+W2B+','+W3B+', \n'
       
    return cbeam_card

    
    







# PBEAM cards =================================================================

# PBEAM,  PID,  MID,     A   ,  I1=Iz  , I2=Iy ,  Iyz   ,  J   ,        ,
#+
#      ,   K1,   K2,         ,         ,       ,        ,      ,        ,
#      ,     ,     ,         ,         ,  N1(A),   N2(A), N1(B),   N2(B),
# PID           -> propery id
# MID           -> material id
# A,I1,I2,I12,J -> Section Properties (I12 default=0)
# K1,K2 -> shear stiffness factor for planes 1 and 2 (0.0 for Bernoulli-Euler bemas)
# N1,N2 -> neutral axis coordinates wrt shear center
# Note: beam local axis is the shear center
# (9 comas per line)
# Same properties all along the beam element


def PBEAM(PID,MID,A,Iy,Iz,J,NA):

    # inputs: 
    # PID,MID   -> integers
    # A,Iy,Iz,J -> reals
    # NA        -> vector (2 components) neutral axis (wrt shear center)


    # Conversion to strings
    PID = str(PID)
    MID = str(MID)
    
    A_A   = str(A)
    I1_A  = str(Iz) #THERE WAS A MISTAKE HERE!
    I2_A  = str(Iy) #THERE WAS A MISTAKE HERE!
    J_A   = str(J)
    K1    =  K2   = '0.0'
    [N1_A, N2_A]  = [str(NA[0]), str(NA[1])]
    [N1_B, N2_B]  = [N1_A, N2_A]
    _ = ' ' #blank/empty

    pbeam_card = ''
    pbeam_card += 'PBEAM'+','+ PID +','+ MID +','+A_A+','+I1_A+','+I2_A+','+  _ +','+ J_A+','+  _ +', \n'\
               +  '+ \n'\
               +  '     '+','+  K1 +','+  K2 +','+ _ +','+ _  +','+ _  +','+  _ +','+  _ +','+  _ +', \n'\
               +  '     '+','+  _  +','+  _  +','+ _ +','+ _  +','+N1_A+','+N2_A+','+N1_B+','+N2_B+', \n'             

    return pbeam_card














# CONM2 cards  ================================================================
# Defines a concetrated mass at a grid point
# CONM2, EID,   G, CID,   M,  X1,  X2, X3,    ,
#      , I11, I21, I22, I31, I32, I33,   ,    ,
# EID      -> element id
# G        -> grid point
# CID      -> coordinate system (default=0)
# M        -> mass
# X1,X2,X3 -> center of gravity of the mass (wrt to G)
# Iij      -> mass moments of inertia wrt center of gravity (1-x, 2-y, 3-z)

def CONM2(EID,G,M,X,Ixx,Iyy,Izz,Ixy,Ixz,Iyz):

    # Inputs
    # EID,G  -> integers
    # M,Iij  -> reals
    # X      -> vector 3 real components
    
    # Conversion to strings
    EID = str(EID)
    G   = str(G)
    M   = str(M)
    [X1,X2,X3]=[str(X[0]),str(X[1]),str(X[2])]
    [I11,I22,I33,I21,I31,I32] = [str(Ixx),str(Iyy),str(Izz),str(Ixy),str(Ixz),str(Iyz)]
    _ = ' ' #blank/empty
           
    # Generation of the card
    conm2_card = ''
    conm2_card += 'CONM2'+','+EID+','+ G +','+ _ +','+ M +','+X1 +','+X2 +','+X3 +','+ _ +', \n'\
               +  '     '+','+I11+','+I21+','+I22+','+I31+','+I32+','+I33+','+ _ +','+ _ +', \n'
      
    return conm2_card






# Eigenvalues request =========================================================
def EIGRL(ND):
# EIGRL, SID, V1, V2, ND 
# SID    -> Set id
# V1,V2  -> Frequency range of interest  (we may include it in the future)
# ND     -> Number of eigenvalues requested
    _ = ' ' #blank/empty
    ND = str(ND)
    eigrl_card = 'EIGRL'+','+'1'+','+_+','+_+','+ND+'\n'
    return eigrl_card









def PROD(PID,MID,A):
    # inputs: 
    # PID,MID   -> integers
    # A         -> reals
    # Conversion to strings
    PID = str(PID)
    MID = str(MID)
    A   = str(A)
    prod_card = ''
    prod_card += 'PROD'+','+ PID +','+ MID +','+ A +'\n'        
    return prod_card


def PSHELL(PID,MID,t):
    # inputs: 
    # PID,MID   -> integers
    # t         -> reals
    # Conversion to strings
    PID = str(PID)
    MID = str(MID)
    t   = str(t)
    pshell_card = ''
    pshell_card += 'PSHELL'+','+ PID +','+ MID +','+ t +'\n'        
    return pshell_card


# RBE2 cards  =================================================================
# EID      -> element id (integer)
# GN       -> independent node (integer)
# CM       -> components (numbers 1,2..6)
# GM       -> dependent nodes (vector of integers)
# Max. 9 items per line 
def RBE2(EID,GN,CM,GM):
    EID = str(EID)
    GN  = str(GN)
    CM  = str(CM)
    rbe2_card = ''
    rbe2_card += 'RBE2'+','+ EID +','+ GN +','+ CM   
    for i,x in enumerate(GM):
        if i != 5:
            rbe2_card += ','+str(x)  
        else:
            rbe2_card += '\n    ,'+str(x)
    rbe2_card += '\n'    
    return rbe2_card



# CORD2R cards  =================================================================
# V_A: position of the origin wrt global SRef
# V_B: position of a point in the z-axis
# V_C: position of a point in the x-axis
def CORD2R(CORDid,V_A,V_B,V_C):
# CORDid      -> id (integer)
# V_A,V_B,V_C -> vectors of 3 components of reals
    _ = ' ' #blank/empty
    CORDid    = str(CORDid)
    [A1,A2,A3]=[str(V_A[0]),str(V_A[1]),str(V_A[2])]
    [B1,B2,B3]=[str(V_B[0]),str(V_B[1]),str(V_B[2])]
    [C1,C2,C3]=[str(V_C[0]),str(V_C[1]),str(V_C[2])]
    cord2r_card  = 'CORD2R'+','+CORDid+','+_+','+A1+','+A2+','+A3+','+B1+','+B2+','+B3+','+'\n'   
    cord2r_card += '      '+','+C1+','+C2+','+C3+'\n'
    return cord2r_card














    
    
# CBAR card ==================================================================
# CBAR, EID, PID,  GA,  GB,  X1,  X2,  X3,    ,
#      ,    ,    , W1A, W2A, W3A, W1B, W2B, W3B,
# EID      -> element id
# PID      -> property id
# GA,GB    -> the two grid points of the beam element
# X1,X2,X3 -> orientation vector of local y-axis (x=local beam axis)
# WA,WB    -> offset vectors from GA,GB of the shear center
def CBAR(EID,PID,GA,GB,X,W):
    
    # inputs: 
    # EID,PID,GA,GB -> integers
    # X,WA,WB -> vectors (3 components)
             
    # Conversion to strings
    EID = str(EID)
    PID = str(PID)
    GA  = str(GA)
    GB  = str(GB)
    [X1,X2,X3]=[str(X[0]),str(X[1]),str(X[2])]
    [W1A,W2A,W3A]=[str(W[0]),str(W[1]),str(W[2])]
    [W1B,W2B,W3B]=[str(W[0]),str(W[1]),str(W[2])]
    _ = ' ' #blank/empty
       
    # Generation of the card
    cbar_card  = ''
    cbar_card += 'CBAR'+','+EID+','+PID+','+GA +','+GB +','+X1 +','+X2 +','+X3 +','+'1'+', \n'\
               +  '   '+','+ _ +','+ _ +','+W1A+','+W2A+','+W3A+','+W1B+','+W2B+','+W3B+', \n'
       
    return cbar_card







# PBAR cards =================================================================

# PBAR,  PID,  MID,     A   ,  I1=Iz  , I2=Iy ,  Iyz   ,  J   ,        ,
#+
#      ,   K1,   K2,         ,         ,       ,        ,      ,        ,
#      ,     ,     ,         ,         ,  N1(A),   N2(A), N1(B),   N2(B),
# PID           -> propery id
# MID           -> material id
# A,I1,I2,I12,J -> Section Properties (I12 default=0)
# K1,K2 -> shear stiffness factor for planes 1 and 2 (0.0 for Bernoulli-Euler bemas)
# N1,N2 -> neutral axis coordinates wrt shear center
# Note: beam local axis is the shear center
# (9 comas per line)
# Same properties all along the beam element
def PBAR(PID,MID,A,Iy,Iz,J):

    # inputs: 
    # PID,MID   -> integers
    # A,Iy,Iz,J -> reals
    # Conversion to strings
    PID = str(PID)
    MID = str(MID)  
    A_A   = str(A)
    I1_A  = str(Iz)
    I2_A  = str(Iy) 
    J_A   = str(J)
    _ = ' ' #blank/empty
    pbar_card = ''
    pbar_card += 'PBAR'+','+ PID +','+ MID +','+A_A+','+I1_A+','+I2_A+','+J_A+','+  _ +', \n'\
           
    return pbar_card




# CRIGD1 cards  =================================================================
# EID      -> element id (integer)
# GN       -> independent node (integer)
# GM       -> dependent nodes (vector of integers)
# Max. 9 items per line 
def CRIGD1(EID,GN,GM):
    EID = str(EID)
    GN  = str(GN)
    crigd1_card = ''
    crigd1_card += 'CRIGD1'+','+ EID +','+ GN
    for i,x in enumerate(GM):
        if i != 6:
            crigd1_card += ','+str(x)  
        else:
            crigd1_card += '\n      ,'+str(x)
    crigd1_card += '\n'    
    return crigd1_card
















