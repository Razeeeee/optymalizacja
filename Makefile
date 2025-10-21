# Makefile dla Windows (projekt Optymalizacja)
# Kompatybilny z MinGW, MSYS2 lub Git Bash

# Zmienne
BUILD_DIR = build
EXECUTABLE = $(BUILD_DIR)/Optymalizacja.exe
CMAKE_BUILD_TYPE ?= Release

# Wykrywanie systemu Windows
ifeq ($(OS),Windows_NT)
    # Windows - używamy cmd
    SHELL = cmd
    RM = del /Q /S
    MKDIR = if not exist "$(BUILD_DIR)" mkdir "$(BUILD_DIR)"
    RM_DIR = if exist "$(BUILD_DIR)" rmdir /S /Q "$(BUILD_DIR)"
    PATH_SEP = \\
    # Sprawdzanie czy plik istnieje w Windows
    FILE_EXISTS = if exist "$(EXECUTABLE)"
    FILE_NOT_EXISTS = if not exist "$(EXECUTABLE)"
else
    # Unix-like (fallback dla MSYS2/Git Bash)
    RM = rm -rf
    MKDIR = mkdir -p $(BUILD_DIR)
    RM_DIR = rm -rf $(BUILD_DIR)
    PATH_SEP = /
    FILE_EXISTS = [ -f "$(EXECUTABLE)" ]
    FILE_NOT_EXISTS = [ ! -f "$(EXECUTABLE)" ]
endif

# Kolory dla czytelniejszego outputu (działają w Git Bash/MSYS2)
GREEN = 
RED = 
YELLOW = 
NC = 

# Sprawdzanie czy terminal wspiera kolory
ifneq ($(TERM),)
    GREEN = \033[0;32m
    RED = \033[0;31m
    YELLOW = \033[1;33m
    NC = \033[0m
endif

.PHONY: all build clean run help dev debug release rebuild install

# Domyślny cel
all: build

# Budowanie projektu
build:
	@echo $(GREEN)Budowanie projektu...$(NC)
	@$(MKDIR)
ifeq ($(OS),Windows_NT)
	@cd $(BUILD_DIR) && cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) ..
	@cd $(BUILD_DIR) && mingw32-make
else
	@cd $(BUILD_DIR) && cmake -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) ..
	@cd $(BUILD_DIR) && make
endif
	@echo $(GREEN)Budowanie zakonczone pomyslnie!$(NC)

# Uruchamianie programu
run:
	@echo $(GREEN)Uruchamianie programu...$(NC)
ifeq ($(OS),Windows_NT)
	@cd $(BUILD_DIR) && Optymalizacja.exe
else
	@cd $(BUILD_DIR) && ./Optymalizacja
endif

# Sprawdzenie czy plik wykonywalny istnieje
check-executable:
ifeq ($(OS),Windows_NT)
	@$(FILE_NOT_EXISTS) ( \
		echo $(YELLOW)Plik wykonywalny nie istnieje. Rozpoczynam budowanie...$(NC) && \
		$(MAKE) build \
	)
else
	@if $(FILE_NOT_EXISTS); then \
		echo "$(YELLOW)Plik wykonywalny nie istnieje. Rozpoczynam budowanie...$(NC)"; \
		$(MAKE) build; \
	fi
endif

# Czyszczenie
clean:
	@echo $(YELLOW)Czyszczenie katalogu build$(PATH_SEP)...$(NC)
	@$(RM_DIR)
	@echo $(GREEN)Czyszczenie zakonczone!$(NC)

# Rebuild (clean + build)
rebuild: clean build

# Debug build
debug:
	@echo $(GREEN)Budowanie w trybie debug...$(NC)
	@$(MKDIR)
ifeq ($(OS),Windows_NT)
	@cd $(BUILD_DIR) && cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Debug ..
	@cd $(BUILD_DIR) && mingw32-make
else
	@cd $(BUILD_DIR) && cmake -DCMAKE_BUILD_TYPE=Debug ..
	@cd $(BUILD_DIR) && make
endif
	@echo $(GREEN)Budowanie debug zakonczone!$(NC)

# Release build
release:
	@echo $(GREEN)Budowanie w trybie release...$(NC)
	@$(MKDIR)
ifeq ($(OS),Windows_NT)
	@cd $(BUILD_DIR) && cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release ..
	@cd $(BUILD_DIR) && mingw32-make
else
	@cd $(BUILD_DIR) && cmake -DCMAKE_BUILD_TYPE=Release ..
	@cd $(BUILD_DIR) && make
endif
	@echo $(GREEN)Budowanie release zakonczone!$(NC)

# Instalacja zależności (dla Windows)
install:
ifeq ($(OS),Windows_NT)
	@echo $(GREEN)Sprawdzanie zainstalowanych narzedzi...$(NC)
	@where cmake >nul 2>&1 || (echo $(RED)BŁĄD: CMake nie jest zainstalowany!$(NC) && exit /b 1)
	@where gcc >nul 2>&1 || where g++ >nul 2>&1 || (echo $(RED)BŁĄD: Kompilator GCC/G++ nie jest zainstalowany!$(NC) && exit /b 1)
	@where mingw32-make >nul 2>&1 || where make >nul 2>&1 || (echo $(RED)BŁĄD: Make nie jest zainstalowany!$(NC) && exit /b 1)
	@echo $(GREEN)Wszystkie wymagane narzedzia sa zainstalowane!$(NC)
else
	@echo $(GREEN)Sprawdzanie zainstalowanych narzedzi...$(NC)
	@command -v cmake >/dev/null 2>&1 || (echo "$(RED)BŁĄD: CMake nie jest zainstalowany!$(NC)" && exit 1)
	@command -v gcc >/dev/null 2>&1 || command -v g++ >/dev/null 2>&1 || (echo "$(RED)BŁĄD: Kompilator GCC/G++ nie jest zainstalowany!$(NC)" && exit 1)
	@command -v make >/dev/null 2>&1 || (echo "$(RED)BŁĄD: Make nie jest zainstalowany!$(NC)" && exit 1)
	@echo "$(GREEN)Wszystkie wymagane narzędzia są zainstalowane!$(NC)"
endif

# Pomoc
help:
	@echo $(GREEN)Dostepne komendy dla Windows:$(NC)
	@echo   $(YELLOW)make build$(NC)     - Buduje projekt w katalogu build$(PATH_SEP)
	@echo   $(YELLOW)make run$(NC)       - Uruchamia program (buduje automatycznie jesli potrzeba)
	@echo   $(YELLOW)make clean$(NC)     - Usuwa katalog build$(PATH_SEP)
	@echo   $(YELLOW)make rebuild$(NC)   - Czysci i buduje od nowa
	@echo   $(YELLOW)make debug$(NC)     - Buduje w trybie debug
	@echo   $(YELLOW)make release$(NC)   - Buduje w trybie release
	@echo   $(YELLOW)make install$(NC)   - Sprawdza czy wszystkie narzedzia sa zainstalowane
	@echo   $(YELLOW)make dev$(NC)       - Buduje i uruchamia program
	@echo   $(YELLOW)make help$(NC)      - Wyswietla te pomoc
	@echo.
	@echo $(GREEN)Wymagania dla Windows:$(NC)
	@echo   - CMake (min. 3.16)
	@echo   - MinGW-w64 lub MSYS2
	@echo   - Git Bash (opcjonalnie, dla kolorow)

# Rozwój (buduje i uruchamia)
dev: build run

# Test - sprawdzenie konfiguracji
test:
	@echo $(GREEN)Test konfiguracji projektu...$(NC)
	@$(MAKE) install
	@$(MAKE) clean
	@$(MAKE) build
	@echo $(GREEN)Test zakończony pomyślnie!$(NC)
