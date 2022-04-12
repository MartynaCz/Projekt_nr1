#!/usr/bin/python
# -*- coding: UTF-8 -*-

import numpy as np
import funkcje as f
import sys
sys.path.append('C:\\Users\\marty\\OneDrive\\Pulpit\\projekt1_inf')

elipsoida_wgs_84 = f.Transformacje(model = "wgs84")

    

#Plik do odczytu
plik = "wsp_inp.txt"
#Odczyt danych z pliku
tablica = np.genfromtxt(plik, delimiter=',', skip_header = 4)

fi_lam_h_N= []
for element in range(0,12):
    fi = elipsoida_wgs_84.XYZ2filamh(tablica[element,0], tablica[element,1], tablica[element,2] )
    fi_lam_h_N.append([fi])
    
    print(fi_lam_h_N)
    
x20_y20 =[]
for element in range(0,12):
     hirvonen = elipsoida_wgs_84.XYZ2filamh(tablica[element,0], tablica[element,1], tablica[element,2] )
     xy= elipsoida_wgs_84.uklad2000(hirvonen[0], hirvonen[1] )
     x20_y20.append(xy)
     print(x20_y20)
     
x92_y92 = []
for element in range(0,12):
     hirvonen = elipsoida_wgs_84.XYZ2filamh(tablica[element,0], tablica[element,1], tablica[element,2] )
     xy= elipsoida_wgs_84.uklad1992(hirvonen[0], hirvonen[1] )
     x92_y92.append(xy)
     print(x92_y92)
     
N_E_U = []     
for element in range(0,12):
    hirvonen = elipsoida_wgs_84.XYZ2filamh(tablica[element,0], tablica[element,1], tablica[element,2] )
    xy= elipsoida_wgs_84.liczNEU(hirvonen[0], hirvonen[1] )
    N_E_U.append(xy)
    print(N_E_U)
    
    



print("Wpisz numer opcji, ktora chcesz wyswietlić\n")
while True:
    print("OPCJE DO WYBORU\n")
    print("0 ---> Wyjście\n")
    print("1 ---> Transformacja współrzędnych kartezjańskich do krzywoliniowych (hirvonen)\n")
    print("2 ---> Transformacja - układ 2000\n")
    print("3 ---> Transformacja - układ 1992\n")
    print("4 ---> Transformacja - współrzedne w trójwymiarze ENU\n")
    print("5 ---> Transformacja współrzędnych krzywoliniowych do kartezjanskich \n")
    print("6 ---> Odleglosc 2D lub 3D\n")
    print("7 ---> Kat azymutu i kat elewacji\n")
    komenda = input()
    if komenda == "0":
        break
    elif komenda == "1":
        fi_lam_h = elipsoida_wgs_84.XYZ2filamh(tablica[0,0],tablica[0,1], tablica[0,2])
    elif komenda == "2":
        x20_y20 = elipsoida_wgs_84.uklad_2000(tablica[0,0], tablica[0,1])
    elif komenda == "3":
        x92_y92 = elipsoida_wgs_84.uklad_1992(tablica[0,0], tablica[0,1])
    elif komenda == "4":
        E_N_U = elipsoida_wgs_84.licz_NEU(tablica[0,0],tablica[0,1], tablica[0,2])
    elif komenda == "5":
        X_Y_Z = elipsoida_wgs_84.filam2XYZ(tablica[0,0],tablica[0,1], tablica[0,2])
   # elif komenda == "6":
   #     funkcje.odl(tablicq)
   # elif komenda == "7":
    #    funkcje.filam2XYZ(tablica)
    else:
        print(" Wybrano opcję, która nie istnieje \n")






#w, r = np.shape(tablica)
#wynik = np.zeros((w, 7))


#i = 0
#for wiersz in tablica:
#    wynik[i,0], wynik[i,1], wynik[i,2] = geo.xyz2plh(wiersz[0], wiersz[1], wiersz[2])
#    wynik[i,3], wynik[i,4] = geo.philam2xy2000(wiersz[i,0], wiersz[i,1])
#    wynik[i,5], wynik[i,6] = geo.philam2xy1992(wiersz[i,0], wiersz[i,1])
    
    
#    i+=1
# zapis: https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.savetxt.html
#np.savetxt("wsp_out.txt", wynik, delimiter=',', fmt = ['%10.2f', '%10.2f', '%10.3f'], header = 'konversja współrzednych geodezyjnych \\ kinga węzka')
