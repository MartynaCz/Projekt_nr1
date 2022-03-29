import math
import numpy as np
from statistics import median
import datetime

def transformacja_kart_do_krzywo(x, y, z):
    '''
    ARGUMENTY:
    X - współrzędna X punktu | typ: float lub int
    Y - współrzędna X punktu | typ: float lub int
    Z - współrzędna X punktu | typ: float lub int

    WYNIKI:
    fi - szerokość geograficzna punktu | typ: float
    lam - długość geograficzna punktu  | typ: float
    h - wysokość punktu               | typ: float
    
    Funkcja, która wykonuje transformacje ze wspolrzednych X, Y, Z do wspolrzednych 
    krzywoliniowych fi, lam, ha oraz zwraca N
    '''
    a = 6378137.000
    e2 = 0.00669438002290
    r = np.sqrt(x**2 + y**2)
    fi_n = math.atan(z / (r*(1-e2)))
    eps = 0.000001/3600*np.pi/180
    fi = fi_n * 2
    while math.sqrt((fi_n - fi)**2) > eps:
        fi = fi_n
        N = a/math.sqrt(1 - e2*math.sin(fi_n)**2)
        h = r/math.cos(fi_n) -N
        fi_n = math.atan(z/(r*(1-e2*(N/(N+h)))))

    N = a / math.sqrt(1 - e2 * math.sin(fi_n) ** 2)
    h = r / math.cos(fi_n) - N
    lam = math.atan(y/x)


    return fi_n, lam, h, N

def zacechowanie_do_00(fi, lam):
    '''
    ARGUMENTY:
    fi- szerokość geograficzna punktu | typ: float
    lam- długość geograficzna punktu  | typ: float
    
    WYNIKI:
    x00- wspolrzedna X | typ: float 
    y00- wspolrzedna X | typ: float    
    Funkcja, która cechuje wspolrzedne krzywoliniowe do ukladu plaskiego 2000
    '''
    a = 6378137
    e2 = 0.00669438002290
    e_2 = e2 / (1 - e2)
    m_0 = 0.999923
    N = a / (math.sqrt(1 - e2 * np.sin(fi) ** 2))
    t = np.tan(fi)
    n2 = e_2 * np.cos(lam) ** 2
    lam = math.degrees(lam)

    if lam > 13.5 and lam < 16.5:
        s = 5
        lam_0 = 15
    elif lam > 16.5 and lam < 19.5:
        s = 6
        lam_0 = 18
    elif lam > 19.5 and lam < 22.5:
        s = 7
        lam_0 = 21
    elif lam > 22.5 and lam < 25.5:
        s = 8
        lam_0 = 24

    lam = math.radians(lam)
    lam_0 = math.radians(lam_0)
    l = lam - lam_0

    A_0 = 1 - (e2 / 4) - (3 * (e2 ** 2)) / 64 - (5 * (e2 ** 3)) / 256
    A_2 = 3 / 8 * (e2 + ((e2 ** 2) / 4) + ((15 * e2 ** 3) / 128))
    A_4 = 15 / 256 * (e2 ** 2 + (3 * (e2 ** 3)) / 4)
    A_6 = (35 * (e2 ** 3)) / 3072

    sigma = a * ((A_0 * fi) - (A_2 * np.sin(2 * fi)) + (A_4 * np.sin(4 * fi)) - (A_6 * np.sin(6 * fi)))

    x = sigma + ((l ** 2) / 2) * (N * np.sin(fi) * np.cos(fi)) * (
                1 + ((l ** 2) / 12) * ((np.cos(fi)) ** 2) * (5 - t ** 2 + 9 * n2 + (4 * n2 ** 2)) + ((l ** 4) / 360) * (
                    (np.cos(fi)) ** 4) * (61 - (58 * (t ** 2)) + (t ** 4) + (270 * n2) - (330 * n2 * (t ** 2))))
    y = l * (N * np.cos(fi)) * (1 + ((((l ** 2) / 6) * (np.cos(fi)) ** 2) * (1 - (t ** 2) + n2)) + (
                ((l ** 4) / (120)) * (np.cos(fi) ** 4)) * (
                                            5 - (18 * (t ** 2)) + (t ** 4) + (14 * n2) - (58 * n2 * (t ** 2))))

    x00 = round(x * m_0, 3)
    y00 = round(y * m_0 + (s * 1000000) + 500000, 3)

    return x00, y00


def srednia(data):
    '''
  
    Funkcja licząca srednia kazdej wspolrzednej: Xsr, Ysr, Zsr | typ: float 
    

    '''
    Sumx = 0
    Sumy = 0
    Sumz = 0
    for sod,X,Y,Z,q in data:
        Sumx += X
        Sumy += Y
        Sumz += Z

    Sx = float(Sumx / len(data))
    Sy = float(Sumy / len(data))
    Sz = float(Sumz / len(data))
    return Sx,Sy,Sz


def licz_NEU(fi, l, X_sr, Y_sr, Z_sr, data):
    '''
    ARGUMENTY:
    fi- szerokość geograficzna punktu | typ: float
    lam- długość geograficzna punktu  | typ: float
    Srednia kazdej wspolrzednej: Xsr, Ysr, Zsr | typ: float     
    WYNIKI:
    Wspolrzedne w trojwymiarowym ukladzie odniesienia U,N,E  | typ: float
    '''
    X=[]
    Y=[]
    Z=[]
    for line in data:
        X.append(line[1])
        Y.append(line[2])
        Z.append(line[3])


    X_sr = float(X_sr)
    Y_sr = float(Y_sr)
    Z_sr = float(Z_sr)

    dX =  []
    dY =  []
    dZ =  []
    NEU = []
    dN =  []
    dE =  []
    dU =  []

    for x, y, z in zip(X, Y, Z):
        delta_X = x - X_sr
        delta_Y = y - Y_sr
        delta_Z = z - Z_sr

        dX.append(delta_X)
        dY.append(delta_Y)
        dZ.append(delta_Z)

    Rt = np.matrix([((-math.sin(fi) * math.cos(l)), (-math.sin(fi) * math.sin(l)), (math.cos(fi))),
                    ((-math.sin(l)), (math.cos(l)), (0)),
                    ((math.cos(fi) * math.cos(l)), (math.cos(fi) * math.sin(l)), (math.sin(fi)))])

    for x, y, z in zip(dX, dY, dZ):
        d = np.matrix([x, y, z])
        d = d.T
        neu = Rt * d
        NEU.append(neu)
        dN.append(float(neu[0]))
        dE.append(float(neu[1]))
        dU.append(float(neu[2]))

    return (NEU, dN, dE, dU)

def regresion(SOD, X):
    '''
    To problem estymacji warunkowej wartości oczekiwanej
    Opisuje współzmienności kilku zmiennych przez dopasowanie do nich funkcji
    
    '''
    A = []
    L = []

    for i in range(len(X)):
        A.append([SOD[i], 1])
        L.append([X[i]])

    A = np.array(A)
    L = np.array(L)

    a = A.T.dot(A)
    a = np.linalg.inv(a)
    a = a.dot(A.T.dot(L))
    residua = a.T.dot(np.vstack([SOD, np.ones(len(SOD))])).T - L
    return a[0], a[1], residua.T[0]

def odchylka(data):
    '''
    Funkcja licząca pierwiastek kwadratowy z wariancji
    
    
    
    '''
    X=[]
    Y=[]
    Z=[]
    for line in data:
        X.append(line[1])
        Y.append(line[2])
        Z.append(line[3])
    sumx = 0
    sumy = 0
    sumz = 0
    Sx, Sy, Sz = srednia(data)
    for i in X:
        sumx += (i-Sx)**2
    for i in Y:
        sumy += (i-Sy)**2
    for i in Z:
        sumz += (i-Sz)**2
    sumx = math.sqrt(sumx / int(len(X)))
    sumy = math.sqrt(sumy / int(len(Y)))
    sumz = math.sqrt(sumz / int(len(Z)))
    return sumx,sumy,sumz

def wariancja(data):
    '''
    Funkcja licząca średnią arytmetyczną kwadratów odchyleń od ich średniej arytmetycznej
    '''
    X=[]
    Y=[]
    Z=[]
    for line in data:
        X.append(line[1])
        Y.append(line[2])
        Z.append(line[3])
    sumx = 0
    sumy = 0
    sumz = 0
    Sx, Sy, Sz = srednia(data)
    for i in X:
        sumx += (i-Sx)**2
    for i in Y:
        sumy += (i-Sy)**2
    for i in Z:
        sumz += (i-Sz)**2
    sumx /= int(len(X))
    sumy /= int(len(Y))
    sumz /= int(len(Z))
    return sumx,sumy,sumz

def odczyt_danych():

    '''
    Funkcja odczyt_danych wczytuje dane z "proj_1_dane.txt" i zwraca
    tabele numpy data(sod,X,Y,Z,Q) oraz dictionary z danymi z nagłówka
    '''

    dictionary = {}
    data = []
    #otwarcie istniejacego pliku
    with open("proj_1_dane.txt") as file:
        for _ in range(8):
            line = file.readline()[2:].replace(":","").replace("\n","").split(maxsplit=1)
            dictionary[line[0]] = line[1]
        for _ in range(2):
            file.readline()
        for line in file:
            if line != "\n":
                line = line.split()
                i = line[1].split(":")
                #wyznaczenie sod
                sod = float(i[0])*3600 + int(i[1])*60 + int(float(i[2]))
                #dodaje wiersze do tabeli
                data.append([sod,float(line[2]),float(line[3]),float(line[4]),float(line[5])])
    data = np.array(data)
    return data, dictionary



def naglowek_raportu(data):

    '''
    Naglowek raportu - funkcja, ktora tworzy nagłówek raportu
    '''

    head = "Imie Nazwisko Autora: Martyna Czarnowska\n" \
           "Data wykonania programu: 27-02-2022\n" \
           f"Data wygenerowania raportu: {datetime.date.today()}\n"
    q = [0,0,0,0,0,0]
    for line in data:
        if line[4] == 1:
            q[0] = q[0] + 1
        if line[4] == 2:
            q[1] = q[1] + 2
        if line[4] == 3:
            q[2] = q[2] + 3
        if line[4] == 4:
            q[3] = q[3] + 4
        if line[4] == 5:
            q[4] = q[4] + 5
        if line[4] == 6:
            q[5] = q[5] + 6
    head += f"Ilość obserwacji z parametrem:\n" \
            f"Q=1 to {int(q[0])}\n" \
            f"Q=2 to {int(q[1]/2)}\n" \
            f"Q=3 to {int(q[2]/3)}\n" \
            f"Q=4 to {int(q[3]/4)}\n" \
            f"Q=5 to {int(q[4]/5)}\n"\
            f"Q=6 to {int(q[5]/6)}\n" \
            f"Łączna suma oberwacji to: {len(data)}\n"
    X, Y, Z = srednia(data)
    X_array = []
    Y_array = []
    Z_array = []
    for line in data:
        X_array.append(line[1])
        Y_array.append(line[2])
        Z_array.append(line[3])
    head += f"Średnia X: {X} | Średnia Y: {Y} | Średnia Z: {Z}\n" \
            f"Mediana X: {median(X_array)} | Mediana Y: {median(Y_array)} | Mediana Z: {median(Z_array)}\n"
    wx,wy,wz = wariancja(data)
    head += f"Wariancja dla X: {wx} | Y: {wy} | Z: {wz} \n"
    head += f"Odchylenie Standardowe dla X: {math.sqrt(wx)} | Y: {math.sqrt(wy)} | Z: {math.sqrt(wz)} \n"
    return head

def stworz_raport(data, what2add):
    '''
    Funkcja stworz_raport służy do generowania raportu
    Generuje dane z wybranych przez użytkownika transformacji oraz dopisuje je do raportu
    W razie nie wybrania przez użytkownika żadnych dodatkowych transformacji generuje raport podstawowy
    '''

    header = naglowek_raportu(data)
    raport_array = data
    raport_array = np.delete(raport_array,4,1)
    titles = "Sod    |       X      |       Y      |       Z      "

    k2k_array = []
    u2k = []
    enu = []
    for i in what2add:
        if i == "2":
            for i in range(len(data)):
                temp = transformacja_kart_do_krzywo(data[i][1], data[i][2], data[i][3])
                k2k_array.append([temp[0],temp[1],temp[2]])
            k2k_array = np.array(k2k_array)
            raport_array = np.append(raport_array, k2k_array, axis=1)
            titles += "|        Fi        |          L         |             H          "
        if i == "3":
            for i in range(len(data)):
                hirvonem = transformacja_kart_do_krzywo(data[i][1], data[i][2], data[i][3])
                temp = zacechowanie_do_00(hirvonem[0], hirvonem[1])
                u2k.append([temp[0],temp[1]])
            u2k = np.array(u2k)
            print(u2k.shape)
            print(raport_array.shape)
            raport_array = np.append(raport_array, u2k, axis=1)
            titles += "|   X_2000   |   Y_2000   "
        if i == "4":
            srednie = srednia(data)
            hirvonem = transformacja_kart_do_krzywo(srednie[0], srednie[1], srednie[2])
            NEU, dN, dE, dU = licz_NEU(hirvonem[0], hirvonem[1], srednie[0], srednie[1], srednie[2], data)
            print("East (E) | North (N) | Up (U)")
            for i in range(len(dN)):
                enu.append([dE[i], dN[i], dU[i]])
            enu = np.array(enu)
            raport_array = np.append(raport_array, enu, axis=1)
            titles += "|      East (E)      |       North (N)      |        Up (U)       "

    raport_array = raport_array.tolist()
    titles += "\n"
    print(header)
    print(titles)
    print(raport_array[0])
    with open("proj_1_raport.txt", 'a') as file:
        file.write(header)
        file.write(titles)
        for i in raport_array:
            file.write(str(i).replace(","," |")[1:-1]+"\n")

def glowna_raport(data):

    '''
    Funkcja glowna_raport - odzielne menu dla generowania raportu
    '''

    what2add = []
    print("Wybierz transformacje która ma zostać dodana do raportu\n"
          "Wybranie ponownie transformacji usunie ją z raportu\n")
    while True:
        print(f"Wybrane transformacje: {what2add}\n"
              f"Opcje:\n"
              f"0 ---> Powrót\n"
              f"1 ---> Generowanie raportu\n"
              f"2 ---> Współrzędne Kartezjańskie do Krzywoliniowych (hirvonen)\n"
              f"3 ---> Transformacja - uklad2000\n"
              f"4 ---> Transformacja do ukladu trojwymiarowego NEU")
        komenda = input()
        if komenda == "0":
            break
        elif komenda == "1":
            if what2add:
                print("Napewno? Wpisz tak lub nie\n")
            else:
                print("Czy chcesz wygenerować podstawowy raport bez innych transformacji? Wpisz tak lub nie\n")
            tln = input().lower()
            if tln == "tak":
                stworz_raport(data, what2add)
                break
            elif tln != "nie":
                print("Zła opcja\n")
        elif komenda == "2":
            if komenda not in what2add:
                what2add.append(komenda)
            else:
                what2add.remove(komenda)
        elif komenda == "3":
            if komenda not in what2add:
                what2add.append(komenda)
            else:
                what2add.remove(komenda)
        elif komenda == "4":
            if komenda not in what2add:
                what2add.append(komenda)
            else:
                what2add.remove(komenda)
