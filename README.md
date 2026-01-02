# gnssdec

Narzędzie do dekodowania sygnałów nawigacji satelitarnej GPS, Galileo i GLONASS z surowych próbek I/Q. Projekt bazuje na GNSS-SDRLIB i został rozszerzony o autorski moduł wyznaczania pozycji, który oblicza współrzędne odbiornika bez potrzeby zewnętrznego oprogramowania. Program obsługuje pełny proces przetwarzania sygnału: akwizycję, śledzenie, dekodowanie wiadomości nawigacyjnych oraz eksport danych w formacie JSON.

## Instalacja zależności na systemie opartym na Debianie:

```sh
sudo apt install build-essential libusb-1.0-0-dev libfec-dev libfftw3-dev
```

## Kompilacja:

```sh
cd bin;
make clean;
make
```

## Sposób uruchomienia:

```sh
./gnssdec -{g/a/l} nagranie.bin
```
