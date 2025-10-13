# Makefile dla projektu Optymalizacja
# Używa CMake do budowania projektu

# Zmienne
BUILD_DIR = build
EXECUTABLE = $(BUILD_DIR)/Optymalizacja
CMAKE_BUILD_TYPE ?= Release

# Kolor dla czytelniejszego outputu
GREEN = \033[0;32m
RED = \033[0;31m
YELLOW = \033[1;33m
NC = \033[0m # No Color

.PHONY: all build clean run help dev

# Domyślny cel
all: build

# Budowanie projektu
build:
	@echo "$(GREEN)Budowanie projektu...$(NC)"
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) ..
	@cd $(BUILD_DIR) && make
	@echo "$(GREEN)Budowanie zakończone pomyślnie!$(NC)"

# Uruchamianie programu (z automatycznym budowaniem jeśli potrzeba)
run: $(EXECUTABLE)
	@echo "$(GREEN)Uruchamianie programu...$(NC)"
	@cd $(BUILD_DIR) && ./Optymalizacja

# Sprawdzenie czy plik wykonywalny istnieje, jeśli nie - buduj
$(EXECUTABLE):
	@if [ ! -f "$(EXECUTABLE)" ]; then \
		echo "$(YELLOW)Plik wykonywalny nie istnieje. Rozpoczynam budowanie...$(NC)"; \
		$(MAKE) build; \
	fi

# Czyszczenie
clean:
	@echo "$(YELLOW)Czyszczenie katalogu build/...$(NC)"
	@rm -rf $(BUILD_DIR)
	@echo "$(GREEN)Czyszczenie zakończone!$(NC)"

# Rebuild (clean + build)
rebuild: clean build

# Debug build
debug:
	@echo "$(GREEN)Budowanie w trybie debug...$(NC)"
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -DCMAKE_BUILD_TYPE=Debug ..
	@cd $(BUILD_DIR) && make
	@echo "$(GREEN)Budowanie debug zakończone!$(NC)"

# Release build (domyślnie)
release:
	@echo "$(GREEN)Budowanie w trybie release...$(NC)"
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -DCMAKE_BUILD_TYPE=Release ..
	@cd $(BUILD_DIR) && make
	@echo "$(GREEN)Budowanie release zakończone!$(NC)"

# Pomoc
help:
	@echo "$(GREEN)Dostępne komendy:$(NC)"
	@echo "  $(YELLOW)make build$(NC)   - Buduje projekt w katalogu build/"
	@echo "  $(YELLOW)make run$(NC)     - Uruchamia program (buduje automatycznie jeśli potrzeba)"
	@echo "  $(YELLOW)make clean$(NC)   - Usuwa katalog build/"
	@echo "  $(YELLOW)make rebuild$(NC) - Czyści i buduje od nowa"
	@echo "  $(YELLOW)make debug$(NC)   - Buduje w trybie debug"
	@echo "  $(YELLOW)make release$(NC) - Buduje w trybie release"
	@echo "  $(YELLOW)make help$(NC)    - Wyświetla tę pomoc"

# Rozwój (buduje i uruchamia)
dev: build run
