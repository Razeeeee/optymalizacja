# Projekt Optymalizacja - Przewodnik budowania

## Czym jest CMake?

**CMake** to narzędzie do budowania projektów C/C++, które automatycznie generuje pliki Makefile lub projekty Visual Studio. Działa na wszystkich systemach (Linux, Windows, macOS) i pozwala używać tego samego kodu konfiguracji wszędzie.

**Wymagania:** CMake, kompilator C++, system budowania (Make/MSBuild)

**Podstawowe komendy:**
```bash
cmake -S . -B build     # Generowanie
cmake --build build     # Budowanie
```

## Czym jest Makefile?

**Makefile** to plik konfiguracyjny dla narzędzia `make`, który definiuje reguły budowania projektu. Zawiera instrukcje, jak kompilować kod źródłowy, linkować biblioteki i zarządzać zależnościami.

### Struktura Makefile:
```makefile
cel: zależności
	komenda_do_wykonania
```

### Podstawowe elementy Makefile:
- **Cele (targets)**: Co chcemy zbudować (np. plik wykonywalny)
- **Zależności**: Od czego zależy dany cel (pliki źródłowe, nagłówki)
- **Komendy**: Jak zbudować dany cel (komendy kompilatora)
- **Zmienne**: Przechowywanie wartości wielokrotnego użytku

### Jak używać Makefile:
```bash
# Uruchomienie domyślnego celu
make

# Uruchomienie konkretnego celu
make clean
make build
make run
```

## Makefile w tym projekcie - najważniejsze aspekty

### 1. **Integracja z CMake**
Nasz Makefile nie kompiluje bezpośrednio plików C++, lecz wykorzystuje CMake jako backend. To zapewnia przenośność między systemami.

### 2. **Zmienne konfiguracyjne**
```makefile
BUILD_DIR = build                    # Katalog budowania
CMAKE_BUILD_TYPE ?= Release          # Domyślny typ budowania
```

### 3. **Phony targets** (cele wirtualne)
```makefile
.PHONY: all build clean run help dev
```
Cele które nie tworzą plików o tych nazwach.

### 4. **Kolorowy output**
Używamy kodów ANSI do kolorowania outputu w terminalu dla lepszej czytelności.

### 5. **Automatyczne zależności**
Cel `run` automatycznie sprawdza czy plik wykonywalny istnieje i buduje go w razie potrzeby.

### 6. **Wsparcie dla różnych trybów budowania**
- `make debug` - kompilacja z informacjami debugowania
- `make release` - kompilacja zoptymalizowana

## Dostępne komendy w projekcie

### Linux/Unix (Makefile)
| Komenda | Opis |
|---------|------|
| `make` lub `make build` | Buduje projekt |
| `make run` | Uruchamia program (buduje automatycznie) |
| `make clean` | Usuwa pliki tymczasowe |
| `make rebuild` | Czyści i buduje od nowa |
| `make debug` | Buduje w trybie debug |
| `make release` | Buduje w trybie release |
| `make dev` | Buduje i uruchamia |
| `make help` | Pokazuje pomoc |

### Windows (Makefile.win)
```bash
# Używanie Makefile dla Windows
make -f Makefile.win build
make -f Makefile.win run
make -f Makefile.win clean

# Lub skopiuj plik i zmień nazwę
copy Makefile.win Makefile
make build
```

| Komenda | Opis |
|---------|------|
| `make -f Makefile.win build` | Buduje projekt używając MinGW |
| `make -f Makefile.win run` | Uruchamia program |
| `make -f Makefile.win install` | Sprawdza wymagane narzędzia |
| `make -f Makefile.win test` | Test pełnej konfiguracji |

## Struktura projektu

```
Optymalizacja/
├── CMakeLists.txt          # Konfiguracja CMake
├── Makefile               # Makefile dla Linux/Unix
├── Makefile.win          # Makefile dla Windows
├── src/                  # Kod źródłowy
│   ├── main.cpp         # Plik główny
│   ├── opt_alg.cpp      # Algorytmy optymalizacji
│   ├── matrix.cpp       # Operacje na macierzach
│   └── ...              # Inne pliki źródłowe
├── build/               # Katalog budowania (generowany)
└── data/                # Dane wejściowe/wyjściowe
```

## Wymagania systemowe

### Linux/Unix:
- **C++17** lub nowszy
- **CMake 3.16** lub nowszy
- **GCC/Clang** (kompilator C++)
- **Make**

### Windows:
- **C++17** lub nowszy
- **CMake 3.16** lub nowszy
- **MinGW-w64** lub **MSYS2** (zawiera GCC i Make)
- **Git Bash** (opcjonalnie, dla kolorowego outputu)

#### Instalacja na Windows:
1. **Pobierz MSYS2** z https://www.msys2.org/
2. **Zainstaluj CMake** z https://cmake.org/download/
3. **W terminalu MSYS2 wykonaj:**
   ```bash
   pacman -S mingw-w64-x86_64-gcc
   pacman -S mingw-w64-x86_64-make
   ```
4. **Dodaj do PATH:** `C:\msys64\mingw64\bin`

---

## 🚀 Tutorial: Jak uruchomić projekt na Windows (krok po kroku)

### Krok 1: Instalacja wymaganych narzędzi

#### Opcja A: MSYS2 (Rekomendowana)
1. **Pobierz i zainstaluj MSYS2:**
   - Idź na https://www.msys2.org/
   - Pobierz installer i uruchom go
   - Zainstaluj w domyślnej lokalizacji `C:\msys64`

2. **Otwórz terminal MSYS2:**
   - Znajdź "MSYS2 MSYS" w menu Start
   - Uruchom terminal

3. **Zainstaluj kompilator i narzędzia:**
   ```bash
   # Aktualizacja pakietów
   pacman -Syu
   
   # Instalacja kompilatora i narzędzi
   pacman -S mingw-w64-x86_64-gcc
   pacman -S mingw-w64-x86_64-make
   pacman -S git
   ```

4. **Pobierz i zainstaluj CMake:**
   - Idź na https://cmake.org/download/
   - Pobierz "Windows x64 Installer"
   - Podczas instalacji **WAŻNE:** zaznacz "Add CMake to system PATH"

#### Opcja B: Git Bash + MinGW (Alternatywna)
1. **Zainstaluj Git for Windows** (zawiera Git Bash)
2. **Zainstaluj MinGW-w64** osobno
3. **Zainstaluj CMake** jak powyżej

### Krok 2: Konfiguracja środowiska

1. **Otwórz terminal MSYS2 MinGW 64-bit** (ważne: nie zwykły MSYS2!)
2. **Sprawdź czy wszystko działa:**
   ```bash
   gcc --version      # Powinno pokazać wersję GCC
   cmake --version    # Powinno pokazać wersję CMake
   make --version     # Powinno pokazać wersję Make
   ```

### Krok 3: Pobranie i budowanie projektu

1. **Sklonuj repozytorium** (jeśli jeszcze nie masz):
   ```bash
   git clone https://github.com/Razeeeee/optymalizacja.git
   cd optymalizacja
   ```

2. **Sprawdź wymagania projektu:**
   ```bash
   make -f Makefile.win install
   ```

3. **Zbuduj projekt:**
   ```bash
   # Pierwsza opcja - używając Makefile dla Windows
   make -f Makefile.win build
   
   # Druga opcja - bezpośrednio CMake
   mkdir build
   cd build
   cmake -G "MinGW Makefiles" ..
   mingw32-make
   ```

4. **Uruchom program:**
   ```bash
   # Opcja 1 - przez Makefile
   make -f Makefile.win run
   
   # Opcja 2 - bezpośrednio
   cd build
   ./Optymalizacja.exe
   ```

### Krok 4: Wygodne użytkowanie

1. **Skopiuj Makefile dla Windows jako domyślny:**
   ```bash
   cp Makefile.win Makefile
   ```

2. **Teraz możesz używać krótkich komend:**
   ```bash
   make build     # Budowanie
   make run       # Uruchamianie
   make clean     # Czyszczenie
   make debug     # Tryb debug
   make help      # Pomoc
   ```

### 🔧 Rozwiązywanie problemów

**Problem:** `cmake: command not found`
- **Rozwiązanie:** Upewnij się, że CMake jest w PATH lub zainstaluj przez MSYS2: `pacman -S mingw-w64-x86_64-cmake`

**Problem:** `gcc: command not found`
- **Rozwiązanie:** Używasz zwykłego terminala MSYS2 zamiast MinGW 64-bit

**Problem:** `make: command not found`
- **Rozwiązanie:** Zainstaluj make: `pacman -S mingw-w64-x86_64-make`

**Problem:** Błędy linkowania
- **Rozwiązanie:** Upewnij się, że używasz generatora "MinGW Makefiles" w CMake

### 💡 Szybkie polecenia (po instalacji)

```bash
# Pełny cykl rozwoju
make clean && make build && make run

# Tylko budowanie i uruchamianie
make dev

# Testowanie różnych trybów
make debug && make run
make release && make run
```

### 🎯 Weryfikacja instalacji

Wykonaj pełny test:
```bash
make -f Makefile.win test
```

Jeśli wszystko przebiegnie pomyślnie, zobaczysz komunikat "Test zakończony pomyślnie!" ✅
