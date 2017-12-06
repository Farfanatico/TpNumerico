import math
NP=97261
CADENCIA=round(-10/10000*(NP-90000)+35)
HC=20
S=0.12224*2*3.141592*12+0.12224**2*2*3.141592
M=7850*3.141592*0.24448*0.01384*(1-0.01384/0.24448)*12
C=480
TINF=round((200/10000)*(NP-90000)+500)+273		

def funcionRadi(T,t):
	ec=-(HC*S*(T-TINF)+float("5.6703e-8")*0.85*S*(T**4-TINF**4))/(M*C)		#derivada
	return ec

def funcionRadiTinfvariable(T,t,Tinf):
	ec=-(HC*S*(T-Tinf)+float("5.6703e-8")*0.85*S*(T**4-Tinf**4))/(M*C)		#derivada
	return ec
	
def funcion(T,t):
	ec=-(HC*S*(T-TINF))/(M*C)		#derivada
	return ec

	
def tp():
	nbol=50
	s=nbol*CADENCIA
	n=50				#s es tiempo total s=nbol*CADENCIA
	t=[]
	T=[]		#temperatura del solido
	for i in range(1,1000):
		T.append(0.0)			#vector nulo
		t.append(0.0) 		
	i=0
	T[i]=293
	t[i]=0
	exportarACSV("\nEuler\n")
	while t[i]<=s:
		T[i+1]=T[i]+(funcion(T[i],t[i]))*CADENCIA    #EULER
		t[i+1]=t[i]+CADENCIA
		exportarACSV(str(t[i]/60)+","+str(T[i]-273))
		i=i+1
	

def RK4():		#cadencia
	nbol=50
	s=nbol*CADENCIA
	t=[]
	T=[]
	for i in range(1,1000):
		T.append(0.0)			#vector nulo
		t.append(0.0) 
	i=0
	T[i]=293
	exportarACSV("\nRK4\n")
	while t[i]<=s:
		q1=CADENCIA*funcion(T[i],t[i])
		q2=CADENCIA*funcion(T[i]+1/2*q1,t[i])
		q3=CADENCIA*funcion(T[i]+1/2*q2,t[i])
		q4=CADENCIA*funcion(T[i]+q3,t[i])
		T[i+1]=T[i]+1/6*(q1+2*q2+2*q3+q4)
		t[i+1]=t[i]+CADENCIA
		exportarACSV(str(t[i]/60)+","+str(T[i]-273))
		i=i+1
	
	
def RK4PUNTO2():
	nbol=50
	s=nbol*CADENCIA
	t=[]
	T=[]
	for i in range(1,1000):
		T.append(0.0)			#vector nulo
		t.append(0.0) 
	i=0
	T[i]=293
	exportarACSVPUNTO2("\nRK4PUNTO2\n")
	while t[i]<=s:
		q1=CADENCIA*funcionRadi(T[i],t[i])
		q2=CADENCIA*funcionRadi(T[i]+1/2*q1,t[i])
		q3=CADENCIA*funcionRadi(T[i]+1/2*q2,t[i])
		q4=CADENCIA*funcionRadi(T[i]+q3,t[i])
		T[i+1]=T[i]+1/6*(q1+2*q2+2*q3+q4)
		t[i+1]=t[i]+CADENCIA
		exportarACSVPUNTO2(str(t[i]/60)+","+str(T[i]-273))
		i=i+1
	aux=0
	c=0
	duracionSK=[]
	exportarACSVPUNTO2("RESULTADOS SOAKING\n")
	for k in range(1,i):
		if (T[i]-T[k]<=10):
			aux+=T[k]-273
			c+=1
			duracionSK.append(t[k]/60)
			exportarACSVPUNTO2(str(t[k]/60)+","+str(T[k]-273))
	exportarACSVPUNTO2("PROMEDIO SOAKING")
	sk=duracionSK[len(duracionSK)-1]-duracionSK[0]
	Tsk=aux/c
	exportarACSVPUNTO2(str(sk)+","+str(Tsk))
	return [Tsk,sk]
	
def analitica():
	tiempo=0
	exportarACSV("\nRESOLUCION ANALITICA\n")
	while (tiempo<=50*CADENCIA):
		T=TINF+(293-TINF)*math.exp(-20*S*tiempo/(M*C))
		exportarACSV (str(tiempo/60)+","+str(T-273))
		tiempo+=CADENCIA
	exportarACSV("\n\n")
			
def punto3():
	nbol=50
	s=nbol*CADENCIA
	t=[]
	T=[]
	for i in range(1,1000):
		T.append(0.0)			#vector nulo
		t.append(0.0) 
	i=0
	T[i]=293
	exportarACSVPUNTO3("\nRK4PUNTO3\n")
	while t[i]<=s:
		if (t[i]<nbol*CADENCIA/2):
			Tvariable=989
		else:
			Tvariable=904
		q1=CADENCIA*funcionRadiTinfvariable(T[i],t[i],Tvariable)
		q2=CADENCIA*funcionRadiTinfvariable(T[i]+1/2*q1,t[i],Tvariable)
		q3=CADENCIA*funcionRadiTinfvariable(T[i]+1/2*q2,t[i],Tvariable)
		q4=CADENCIA*funcionRadiTinfvariable(T[i]+q3,t[i],Tvariable)
		T[i+1]=T[i]+1/6*(q1+2*q2+2*q3+q4)
		t[i+1]=t[i]+CADENCIA
		exportarACSVPUNTO3(str(t[i]/60)+","+str(T[i]-273))
		i=i+1
	
	exportarACSVPUNTO3("RESULTADOS SOAKING\n")
	for k in range(1,i):
		if (T[i]-T[k]<10):
			exportarACSVPUNTO3(str(t[k]/60)+","+str(T[k]-273))
	
def RK4PUNTO4(T1semilla,T2semilla,exportartemp=False):
	nbol=50
	s=nbol*CADENCIA
	t=[]
	T=[]
	for i in range(1,1000):
		T.append(0.0)			#vector nulo
		t.append(0.0) 
	i=0
	T[i]=293
	while t[i]<=s:
		if (t[i]<s/2):
			Tvariable=T1semilla
		else:
			Tvariable=T2semilla
		q1=CADENCIA*funcionRadiTinfvariable(T[i],t[i],Tvariable)
		q2=CADENCIA*funcionRadiTinfvariable(T[i]+1/2*q1,t[i],Tvariable)
		q3=CADENCIA*funcionRadiTinfvariable(T[i]+1/2*q2,t[i],Tvariable)
		q4=CADENCIA*funcionRadiTinfvariable(T[i]+q3,t[i],Tvariable)
		T[i+1]=T[i]+1/6*(q1+2*q2+2*q3+q4)
		t[i+1]=t[i]+CADENCIA
		if (exportartemp):
			exportarACSVPUNTO4(str(t[i]/60)+","+str(T[i]-273))
		i=i+1
	aux=0
	c=0
	duracionSK=[]
	for k in range(0,i):
		if (T[i]-T[k]<=10):
			aux+=T[k]-273
			c+=1
			duracionSK.append(t[k]/60)
	sk=duracionSK[len(duracionSK)-1]-duracionSK[0]
	Tsk=aux/c
	if (exportartemp):
		for k in range(1,i):
			if (T[i]-T[k]<10):
				exportarACSVPUNTO4(str(t[k]/60)+","+str(T[k]-273))
		exportarACSVPUNTO4("\nPROMEDIO TIEMPO SK\n")
		exportarACSVPUNTO4(str(Tsk)+","+str(sk)+"\n")
	return [Tsk,sk]	
	

		
def punto4(Tskobj):
	skobj=10
	Tz=[Tskobj,Tskobj]
	Tz  =[Tz[0],Tz[1]]
	for i in range(1,100):
		[Tsk,sk] = RK4PUNTO4(Tz[0]+273,Tz[1]+273)
		F2 = sk
		F1= Tsk                  
		TZ_new =[ Tz[0]-0.25*(F1-Tskobj)-0.75*(F2-skobj),Tz[1]-0.75*(F1-Tskobj)-0.25*(F2-skobj)] #                       
		exportarACSVPUNTO4 (str(TZ_new[0])+","+str(TZ_new[1])+","+str(i))		
		if (abs(TZ_new[0]-Tz[0])/TZ_new[0] < 0.5e-3 and abs((TZ_new[1]-Tz[1])/TZ_new[1])< 0.5e-3):
			Tz = TZ_new
			break	
		Tz = TZ_new
	exportarACSVPUNTO4("\n\n RESULTADOS TEMPERATURA PARA SOAKING OBTENIDO\n")
	RK4PUNTO4(Tz[0]+273,Tz[1]+273,True)

def exportarACSV(texto):
	archivo=open('resultados.csv','a')
	archivo.write(str(texto+'\n'))
	archivo.close()

def exportarACSVPUNTO2(texto):
	archivo=open('resultadospunto2.csv','a')
	archivo.write(str(texto+'\n'))
	archivo.close()

def exportarACSVPUNTO3(texto):
	archivo=open('resultadospunto3.csv','a')
	archivo.write(str(texto+'\n'))
	archivo.close()
	
def exportarACSVPUNTO4(texto):
	archivo=open('resultadospunto4.csv','a')
	archivo.write(str(texto+'\n'))
	archivo.close()
	
def limpiarcsv():
	archivo=open('resultadospunto4.csv','w')
	archivo.close()
	archivo=open('resultadospunto3.csv','w')
	archivo.close()
	archivo=open('resultados.csv','w')
	archivo.close()
	archivo=open('resultadospunto2.csv','w')
	archivo.close()
		
def main():	
	limpiarcsv()	
	tp()		
	RK4()
	analitica()
	RK4PUNTO2()
	punto3()
	punto4(626.02)
	punto4(622.61)
	punto4(672.61)

main()
