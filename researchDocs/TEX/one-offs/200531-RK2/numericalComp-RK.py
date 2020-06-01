"""
File meant to show numerical integration methods applied via python
Structured in a way that is more related to the simulation method in PSLTDSim

lambda is the python equivalent of matlab anonymous functions

5/31/20 Addition of RK2 (heun) and thrid order system

"""
# Package Imports
import numpy as np
import matplotlib.pyplot as plt

# Method Definitions
def euler(fp, x0, y0, ts):
    """
    fp = Some derivative function of x and y
    x0 = Current x value
    y0 = Current y value
    ts = time step
    Returns y1 using Forward Euler or tangent line method
    """    
    return y0 + fp(x0,y0)*ts

def adams2(fp, x0, y0, xN, yN, ts):
    """
    fp = Some derivative function of x and y
    x0 = Current x value
    y0 = Current y value    
    xN = Previous x value
    yN = Previous y value
    ts = time step
    Returns y1 using Adams-Bashforth two step method
    """
    return y0 + (1.5*fp(x0,y0) - 0.5*fp(xN,yN))*ts 

def rk4(fp, x0, y0, ts):
    """
    fp = Some derivative function of x and y
    x0 = Current x value
    y0 = Current y value
    ts = time step
    Returns y1 using Runge-Kutta method
    """
    k1 = fp(x0, y0)
    k2 = fp(x0 +ts/2, y0+ts/2*k1)
    k3 = fp(x0 +ts/2, y0+ts/2*k2)
    k4 = fp(x0 +ts, y0+ts*k3)    
    return y0 + ts/6*(k1+2*k2+2*k3+k4)

def rk2(fp, x0, y0, ts):
    """
    fp = Some derivative function of x and y
    x0 = Current x value
    y0 = Current y value
    ts = time step
    Returns y1 using Runge-Kutta 2 method
    """
    k1 = fp(x0, y0)
    k2 = fp(x0+ts, y0+k1*ts) 
    return y0 + ts/2*(k1+k2)

def heun(fp, x0, y0, ts):
    """
    fp = Some derivative function of x and y
    x0 = Current x value
    y0 = Current y value
    ts = time step
    Returns y1 using Heun's method
    """
    fp1 = fp(x0, y0)
    y1E = y0 + fp1*ts
    pf2 = fp(x0+ts, y1E) 
    return y0 + ts/2*(fp1 + pf2)

def trapezoidalPost(x,y):
    """
    x = list of x values
    y = list of y values
    Returns integral of y over x.
    Assumes full lists / ran post simulation
    """
    integral = 0
    for ndx in range(1,len(x)):
        integral+= (y[ndx]+y[ndx-1])/2 * (x[ndx]-x[ndx-1])
    return integral
    

# Case Selection
for caseN in range(0,4):
    blkFlag = False # for holding plots open
    if caseN == 0:
        # Trig example
        caseName = 'Sinusoidal Example'
        tStart =0
        tEnd = 3
        numPoints = 6*2

        ic = [0,0] # initial condition x,y
        fp = lambda x, y: -2*np.pi*np.cos(2*np.pi*x)
        f = lambda x,c: -np.sin(2*np.pi*x)+c
        findC = lambda x,y: y+np.sin(2*np.pi*x)
        c = findC(ic[0],ic[1])
        calcInt = ( 1/(2*np.pi)*np.cos(2*np.pi*tEnd)+c*tEnd -
                    1/(2*np.pi)*np.cos(2*np.pi*ic[0])-c*ic[0] ) 

    elif caseN == 1:
        # Exp example
        caseName = 'Exponential Example'
        tStart =0
        tEnd = 3
        numPoints = 3

        ic = [0,0] # initial condition x,y
        fp = lambda x, y: np.exp(x)
        f = lambda x,c: np.exp(x)+c
        findC = lambda x, y: y-np.exp(x)
        c= findC(ic[0],ic[1])
        calcInt = np.exp(tEnd)+c*tEnd-np.exp(ic[0])+c*ic[0]

    elif caseN == 2:
        # Log example
        caseName = 'Logarithmic Example'
        tStart =1
        tEnd = 4
        numPoints = 3
        blkFlag = False # for holding plots open

        ic = [1,1] # initial condition x,y
        fp = lambda x, y: 1/x
        f = lambda x,c: np.log(x)+c
        findC = lambda x, y: y-np.log(x)
        c= findC(ic[0],ic[1])
        calcInt = (tEnd*np.log(tEnd)- tEnd +c*tEnd -
                   ic[0]*np.log(ic[0])+ ic[0] -c*ic[0])
    else:
        # step multi order system
        caseName = 'Step Input Third Order System Example'
        tStart =0
        tEnd = 5
        numPoints = 5*2
        blkFlag = True # for holding plots open

        U = 1
        T0 = 0.4
        T2 = 4.5
        T1 = 5
        T3 = -1
        T4 = 0.5

        alphaNum = (T1*T3)
        alphaDen = (T0*T2*T4)
        alpha = alphaNum/alphaDen

        num = alphaNum*np.array([1, 1/T1+1/T3, 1/(T1*T3)])
        den  = alphaDen*np.array([1, 1/T4+1/T0+1/T2, 1/(T0*T4)+1/(T2*T4)+1/(T0*T2),
                         1/(T0*T2*T4)])
    
        # PFE
        A = ((1/T1-1/T0)*(1/T3-1/T0))/((1/T2-1/T0)*(1/T4-1/T0))
        B = ((1/T1-1/T2)*(1/T3-1/T2))/((1/T0-1/T2)*(1/T4-1/T2))
        C = ((1/T1-1/T4)*(1/T3-1/T4))/((1/T0-1/T4)*(1/T2-1/T4))
        
        initState = 0 # for steady state start
        ic = [0,0] # initial condition x,y
        fp = lambda x, y: alpha*(A*np.exp(-x/T0)+B*np.exp(-x/T2)+C*np.exp(-x/T4))
        f = lambda x, c: alpha*(-T0*A*np.exp(-x/T0)-T2*B*np.exp(-x/T2)-T4*C*np.exp(-x/T4))+c
        findC = lambda x, y : alpha*(A*T0+B*T2+C*T4)
        c = findC(ic[0], ic[1])
        calcInt = (
            alpha*A*T0**2*np.exp(-tEnd/T0) +
            alpha*B*T2**2*np.exp(-tEnd/T2) +
            alpha*C*T4**2*np.exp(-tEnd/T4) +
            c*tEnd -
            alpha*(A*T0**2+B*T2**2+C*T4**2)
            )# Calculated integral

    # Initialize current value dictionary
    # Shown to mimic PSLTDSim record keeping
    cv={
        't' :ic[0],
        'yE': ic[1],
        'yH' : ic[1],
        'yRK2': ic[1],
        'yRK4': ic[1],    
        'yAB': ic[1],
        }

    # Calculate time step
    ts = (tEnd-tStart)/numPoints
    # Find C from integrated equation for exact soln
    c = findC(ic[0], ic[1])
    # Calculate exact solution
    tExact = np.linspace(tStart,tEnd, 10000)
    yExact = f(tExact, c)

    # Initialize running value lists
    t=[]
    yE=[]
    yH = []
    yRK2 =[]
    yRK4 =[]
    yAB = []

    t.append(cv['t'])
    yE.append(cv['yE'])
    yH.append(cv['yH'])
    yRK2.append(cv['yRK2'])
    yRK4.append(cv['yRK4'])
    yAB.append(cv['yAB'])  

    # Start Simulation
    while cv['t']< tEnd:

        # Calculate Euler result
        cv['yE'] = euler( fp, cv['t'], cv['yE'],  ts )
        # Calculate Runge-Kutta results
        cv['yRK2'] = rk2( fp, cv['t'], cv['yRK2'],  ts )
        cv['yRK4'] = rk4( fp, cv['t'], cv['yRK4'],  ts )
        # Calculate Heun's results
        cv['yH'] = heun( fp, cv['t'], cv['yH'],  ts )

        # Calculate Adams-Bashforth result
        if len(t)>=2:
            cv['yAB'] = adams2( fp, cv['t'], cv['yAB'], t[-2], yAB[-2], ts )
        else:
            # Required to handle first step when a -2 index doesn't exist
            cv['yAB'] = adams2( fp, cv['t'], cv['yAB'], t[-1], yAB[-1], ts )

            
        # Log calculated results
        yE.append(cv['yE'])
        yH.append(cv['yH'])
        yRK2.append(cv['yRK2'])
        yRK4.append(cv['yRK4'])
        yAB.append(cv['yAB'])
    
        # Increment and log time
        cv['t'] += ts
        t.append(cv['t'])

    # Generate Plot
    fig, ax = plt.subplots()
    ax.set_title('Approximation Comparison\n' + caseName)
    
    #Plot all lines
    ax.plot(tExact,yExact,
            c=[0,0,0],
            linewidth=2,
            label="Exact")
    ax.plot(t,yE,
            marker='o',
            fillstyle='none',
            linestyle=':',
            c=[0.7,0.7,0.7],
            label="Euler")
    ax.plot(t,yRK4,
            marker='*',
            markersize=10,
            fillstyle='none',
            linestyle=':',
            c=[1,0,1],
            label="RK4")
    ax.plot(t,yRK2,
            marker='x',
            markersize=8,
            fillstyle='none',
            linestyle=':',
            c=[0,0,1],
            label="RK2")
    ax.plot(t,yH,
            marker='+',
            markersize=8,
            fillstyle='none',
            linestyle=':',
            c=[1, .5, 0],
            label="Heun")
    ax.plot(t,yAB,
            marker='s',
            fillstyle='none',
            linestyle=':',
            c =[0,1,0],
            label="AB2")

    # Format Plot
    fig.set_dpi(150)
    fig.set_size_inches(9, 2.5)
    ax.set_xlim(min(t), max(t))
    ax.grid(True, alpha=0.25)
    ax.legend(loc='best',  ncol=3)
    ax.set_ylabel('y Value')
    ax.set_xlabel('x Value')
    fig.tight_layout()    
    plt.show(block = blkFlag)
    plt.pause(0.00001)

    # Trapezoidal Integration
    exactI = trapezoidalPost(tExact,yExact)
    Eint = trapezoidalPost(t,yE)
    Hint = trapezoidalPost(t,yH)
    RKint = trapezoidalPost(t,yRK4)
    RK2int = trapezoidalPost(t,yRK2)
    ABint = trapezoidalPost(t,yAB)

    print("\n%s" % caseName)
    print("time step:  %.2f" % ts)
    print("Method: Trapezoidal Int\t Absolute Error from calculated")
    print("Calc: \t%.9f\t%.9f" % (calcInt ,abs(calcInt-calcInt)))
    print("Exact: \t%.9f\t%.9f" % (exactI ,abs(calcInt-exactI)))
    print("RK4: \t%.9f\t%.9f" % (RKint,abs(calcInt-RKint)))
    print("RK2: \t%.9f\t%.9f" % (RK2int,abs(calcInt-RK2int)))
    print("Heun: \t%.9f\t%.9f" % (RK2int,abs(calcInt-Hint)))
    print("AB2: \t%.9f\t%.9f" % (ABint,abs(calcInt-ABint)))
    print("Euler: \t%.9f\t%.9f" % (Eint,abs(calcInt-Eint)))




