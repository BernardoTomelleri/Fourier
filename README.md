# Fusofourier
Esercizio obbligatorio sull'espansione in serie di Fourier:

Per la compilazione la memoria base di LaTeX potrebbe non essere sufficiente:\
Per gli utenti di TexLive è sufficiente modificare il file di configurazione texmf.cnf (di default: C:\texlive\2019\texmf.cnf)
individuabile nel vostro sistema da riga di comando tramite `kpsewhich -a texmf.cnf`.\
Aggiungendo le righe:
```
main_memory = 10000000
save_size  = 10000000
```
(al di sotto di)
```
OSFONTDIR = $SystemRoot/fonts//
```
(Valori totalmente arbitrari, di default 5000000 e 100000 rispettivamente)\
Questi valori saranno usati dal compilatore una volta ricostruiti i file di formato,
da terminale con privilegi amministrativi, mediante il comando: `fmtutil-sys --all`\
Si veda https://tex.stackexchange.com/questions/7953/how-to-expand-texs-main-memory-size-pgfplots-memory-overload
per ulteriori informazioni.
