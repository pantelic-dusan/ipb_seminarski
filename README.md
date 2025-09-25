# Istraživački projekat
### Tema: Utvrdjivanje potencijalne veze između preslikavanja kodona u amino kiselinu u uređenom i neuređenom delu proteina
### Predmet: Istraživanje podataka u bioinformatici
### Autor: Dušan Pantelić 1062/2024

## Opis projekta
Cilj ovog projekta je da se istraži potencijalna veza između preslikavanja kodona u amino kiselinu u uređenom i neuređenom delu proteina.
Projekat obuhvata prikupljanje podataka iz relevantnih baza podataka, njihovu obradu i analizu kako bi se identifikovale moguće korelacije između strukturalnih karakteristika proteina i njihovih genetskih kodova.

#### Korišćeni podaci
- ***Escherichia coli str. K-12 substr. MG1655, complete genome***(https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3)

#### Korišćeni alati:
- ***IUpred3*** (https://iupred3.elte.hu/)
- ***IsUnstruct*** (http://bioinfo.protres.ru/IsUnstruct/)

#### Korišćene biblioteke:
- ***Biopython*** (https://biopython.org/)

## Uputstvo za pokretanje
Napomena: Koriscen je Python 3.12. Preporučuje se kreiranje virtuelnog okruženja.

### 1. Kloniranje repozitorijuma
```bash 
  git clone https://github.com/pantelic-dusan/ipb_seminarski.git
```

### 2.  Instalacija potrebnih biblioteka
```bash 
  cd ./ipb_seminarski
  pip install -r requirements.txt
```

### 3. Prikupljanje i obrada podataka
Koristiti generisan fajl `data/final_cds_data.json` ili generisati novi.
Za generisanje fajla otvoriti `Prikupljanje_i_Obrada_Podataka.ipynb` u Jupyter Notebook-u i pokrenuti sve ćelije.
