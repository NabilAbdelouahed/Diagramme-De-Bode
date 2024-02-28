import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as nppol

def gain(complexe):
    g=np.abs(complexe)
    g=20*np.log10(g)
    return(g)
    
def phase(complexe):
    return(np.angle(complexe)*(180/(np.pi)))
    
def bodegeneral(num,denum):
    if denum != [0] :
        
        #préparation des valeurs de w
        
        tabnum=np.array(num)
        tabdenum=np.array(denum)
        p10=np.linspace(-4,4,10000)
        w=np.power(10,p10)
        
        #déterminer les gains et les phases
        
        k=[]
        n=[]
        h=[] 
        g=[] 
        pha=[]
        for i in range(len(w)):
            k.append(nppol.polyval(w[i]*1j,tabnum))
            n.append(nppol.polyval(w[i]*1j,tabdenum))
            h.append(k[i]/n[i])
            g.append(gain(h[i]))
            pha.append(phase(h[i]))
            
        #traçage
        
        font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,}
        plt.figure()
        
        #courbe du gain
        
        plt.subplot(211)
        plt.plot(w,g,  linewidth=3  )
        plt.semilogx()
        plt.grid(True)
        plt.xlabel("pulsation (rad/s)" , fontdict=font) 
        plt.ylabel("gain (db)" , fontdict=font)
        
        #courbe de la phase

        plt.subplot(212)
        plt.plot(w,pha, linewidth=3)
        plt.semilogx()
        plt.grid(True)
        plt.show()
        plt.xlabel("pulsation (rad/s)" , fontdict=font) 
        plt.ylabel("déphasage (degrés)" , fontdict=font)
        
        
    return(None)

def bodepremierordre(k,tau,denum):
    
    #distinguer les deux cas possibles
    
    def ft(p):
        if denum==False:
            f=k*(1+tau*p)
        else :
            f=k/(1+tau*p)
        return(f)
    
    #préparation des valeurs de w et calcul des phases et gains
    
    p=np.linspace(-4,4,100000)
    w=np.power(10,p)
    h=ft(w*1j)
    G=gain(h)
    P=phase(h)
    
    #traçage
    
    font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,}
    plt.figure()
    
    #courbe du gain
    
    plt.subplot(211)
    plt.plot(w,G, label="courbe" , linewidth=3  )
    plt.semilogx()
    plt.grid(True)
    
    plt.xlabel("pulsation (rad/s)" , fontdict=font) 
    plt.ylabel("gain (db)" , fontdict=font)
    
    #traçage des assymptotes du gain
    
    plt.plot([0 , (1/tau)],[(gain(k)) , (gain(k))],"r--", label="tracé asymptotique" , linewidth=3 )
    
    plt.plot([(1/tau), 10**4],[(gain(k)) , (gain(h[-1]))],"r--" , linewidth=3 )
    
    if denum==False:
        plt.text(10**3,gain(ft(10**3*1j))-20,'+20db/dec' , fontsize=14, color='g')
    else:
        plt.text(10**3,gain(ft(10**3*1j))+20,'-20db/dec' , fontsize=14, color='g')
    
    plt.legend()
    
    #traçage de la phase
    
    plt.subplot(212)
    plt.plot(w,P, label="courbe", linewidth=3) 
    plt.semilogx()
    plt.grid(True)
    
    plt.xlabel("pulsation (rad/s)" , fontdict=font) 
    plt.ylabel("déphasage (degrés)" , fontdict=font)
    
    #traçage des assymptotes de la phase
    
    plt.plot([0, (1/tau)],[0 , 0],"r--" ,label="tracé asymptotique", linewidth=3 )
    
    if denum==False:
        
        plt.plot([(1/tau), (1/tau)],[0 , 90],"r--" , linewidth=3 )
        
        plt.plot([(1/tau), 10**4],[90 , 90],"r--" , linewidth=3 )
    else:
        
        plt.plot([(1/tau), (1/tau)],[0 , -90],"r--" , linewidth=3 )
        
        plt.plot([(1/tau), 10**4],[-90 , -90],"r--" , linewidth=3 )
        
    
    plt.legend()
    plt.show()
    
    return(None)

def bodesecondordre(k,ksi,w0,denum):
    
    #distinguer les deux cas possibles
    
    def ft(p):
        if denum==False:
            f=k*(1+(2*ksi/w0)*p+(1/w0**2)*p**2)
        else :
            f=k/(1+(2*ksi/w0)*p+(1/w0**2)*p**2)
        return(f)
    
    #préparation des valeurs de w et calcul des phases et gains
    
    p=np.linspace(-4,4,100000)
    w=np.power(10,p)
    h=ft(w*1j)
    G=gain(h)
    P=phase(h)
    
    #traçage
    
    font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,}
    plt.figure()
    
    #courbe du gain
    
    plt.subplot(211)
    plt.plot(w,G, label="courbe" , linewidth=3  )
    plt.semilogx()
    plt.grid(True)
    
    plt.xlabel("pulsation (rad/s)" , fontdict=font) 
    plt.ylabel("gain (db)" , fontdict=font)
    
    #traçage des assymptotes du gain
    
    if ksi<=1:
            
        plt.plot([0 , w0],[(gain(k)) , (gain(k))],"r--", label="tracé asymptotique" , linewidth=3 )
        
        plt.plot([w0 , 10**4],[(gain(k)) , G[-1]],"r--", linewidth=3 )
        
        if denum==False:
            plt.text(10**3,gain(ft(10**3*1j))-30,'+40db/dec' , fontsize=14, color='g')
        else:
            plt.text(10**3,gain(ft(10**3*1j))+30,'-40db/dec' , fontsize=14, color='g')
    
    
    elif ksi>1:
        
        t1=(w0/(ksi+((ksi**2)-1)**(1/2)))  
        t2=(w0/(ksi-((ksi**2)-1)**(1/2)))
            
        if denum==False: 
            
            plt.plot([0 , t1],[(gain(k)) , (gain(k))],"r--", label="tracé asymptotique" , linewidth=3 )
            
            plt.plot([t1 , t2],[(gain(k)) , gain(k)+gain(t2)-gain(t1)],"r--", linewidth=3 )
            
            plt.plot([t2 , 10**4],[gain(k)+gain(t2)-gain(t1) , G[-1]],"r--", linewidth=3 )
            
            plt.plot([t1 , w0],[(gain(k)) , gain(k)], "c--" , linewidth=3 )
            
            plt.plot([w0 , t2],[(gain(k)) , gain(k)+gain(t2)-gain(t1)], "c--" , linewidth=3 )
            
            plt.text(w0/ksi,gain(ft(w0/ksi*1j))+20,'+20db/dec' , fontsize=14, color='g')
            plt.text(10**3,gain(ft(10**3*1j))+30,'+40db/dec' , fontsize=14, color='g')
        
        else:
            
            plt.plot([0 , t1],[(gain(k)) , (gain(k))],"r--", label="tracé asymptotique" , linewidth=3 )
            
            plt.plot([t1 , t2],[(gain(k)) , gain(k)-gain(t2)+gain(t1)],"r--", linewidth=3 )
            
            plt.plot([t2 , 10**4],[gain(k)-gain(t2)+gain(t1) , G[-1]],"r--", linewidth=3 )
            
            plt.plot([t1 , w0],[(gain(k)) , gain(k)], "c--" , linewidth=3 )
            
            plt.plot([w0 , t2],[gain(k) , gain(k)-gain(t2)+gain(t1)], "c--" , linewidth=3 )
            
            plt.text(w0/ksi,gain(ft(w0/ksi*1j))-20,'-20db/dec' , fontsize=14, color='g')
            plt.text(10**3,gain(ft(10**3*1j))-30,'-40db/dec' , fontsize=14, color='g')
            
    plt.legend()
    
    #traçage de la phase
    
    plt.subplot(212)
    plt.plot(w,P, label="courbe", linewidth=3) 
    plt.semilogx()
    plt.grid(True)
    
    plt.xlabel("pulsation (rad/s)" , fontdict=font) 
    plt.ylabel("déphasage (degrés)" , fontdict=font)
    
    #traçage des asymptotes du déphasage
    
    if denum==False :
        if ksi<=1/(2**(1/2)):
            
            plt.plot([0 , w0],[0 , 0],"r--", label="tracé asymptotique" , linewidth=3 )
            
            plt.plot([w0 , w0],[0 , 180],"r--", linewidth=3 )
            
            plt.plot([w0 , 10**4],[180 , 180],"r--", linewidth=3 )
        
        elif ksi>1:
            
            t1=(w0/(ksi+((ksi**2)-1)**(1/2)))  
            t2=(w0/(ksi-((ksi**2)-1)**(1/2)))
            
            plt.plot([0 , t1],[0 , 0],"r--", label="tracé asymptotique" , linewidth=3 )
            
            plt.plot([t1 , t1],[0 , 90],"r--", linewidth=3 )
            
            plt.plot([t1 , t2],[90 , 90],"r--", linewidth=3 )
            
            plt.plot([t2 , t2],[90 , 180],"r--" , linewidth=3 )
            
            plt.plot([t2 , 10**4],[180 , 180],"r--", linewidth=3 )
            
    else:
        if ksi<=1/(2**(1/2)):
            
            plt.plot([0 , w0],[0 , 0],"r--", label="tracé asymptotique" , linewidth=3 )
            
            plt.plot([w0 , w0],[0 , -180],"r--", linewidth=3 )
            
            plt.plot([w0 , 10**4],[-180 , -180],"r--", linewidth=3 )
       
        elif ksi>1:
            
            t1=(w0/(ksi+((ksi**2)-1)**(1/2)))  
            t2=(w0/(ksi-((ksi**2)-1)**(1/2)))
            
            plt.plot([0 , t1],[0 , 0],"r--", label="tracé asymptotique" , linewidth=3 )
            
            plt.plot([t1 , t1],[0 , -90],"r--", linewidth=3 )
            
            plt.plot([t1 , t2],[-90 , -90],"r--", linewidth=3 )
            
            plt.plot([t2 , t2],[-90 , -180],"r--" , linewidth=3 )
            
            plt.plot([t2 , 10**4],[-180 , -180],"r--", linewidth=3 )
            
            
    plt.legend()
    plt.show()
    
    return(None)
bodesecondordre(2,5,12,True)