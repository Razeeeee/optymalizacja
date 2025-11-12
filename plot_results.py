import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# Dane z tabeli
data = {
    'a': [4, 4.4934, 5],
    'zewnetrzna': {
        'x1': [2.2061079, 2.0261231, 1.9985517],
        'x2': [3.291583, 3.9436059, 3.9467881],
        'r': [4.0243699, 4.4757189, 4.4913443],
        'y': [-0.19191422, -0.21713976, -0.21678323],
        'calls': [53.44, 30.76, 22.02]
    },
    'wewnetrzna': {
        'x1': [2.1750799, 2.3838197, 2.609011],
        'x2': [2.1581446, 2.3990589, 2.6106106],
        'r': [3.0650645, 3.3825995, 3.690866],
        'y': [0.0262363063, -0.069397237, -0.1414027],
        'calls': [61.74, 62.51, 58.22]
    }
}

# Funkcja celu: f(x1, x2) = sin(pi*r)/(pi*r), gdzie r = sqrt((x1/pi)^2 + (x2/pi)^2)
def objective_function(x1, x2):
    r_arg = np.sqrt((x1 / np.pi)**2 + (x2 / np.pi)**2)
    # Unikamy dzielenia przez zero
    f = np.where(r_arg < 1e-10, 1.0, np.sin(np.pi * r_arg) / (np.pi * r_arg))
    return f

# Utworzenie siatki dla wykresu poziomic
x1_range = np.linspace(0, 6, 400)
x2_range = np.linspace(0, 6, 400)
X1, X2 = np.meshgrid(x1_range, x2_range)
Z = objective_function(X1, X2)

# Tworzenie wykresów dla każdego parametru a
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
fig.suptitle('Optymalizacja z funkcjami kary - wyniki dla różnych wartości parametru a', fontsize=14, fontweight='bold')

for idx, a_val in enumerate(data['a']):
    ax = axes[idx]
    
    # Wykres poziomic funkcji celu
    contour = ax.contour(X1, X2, Z, levels=20, cmap='viridis', alpha=0.6)
    contourf = ax.contourf(X1, X2, Z, levels=20, cmap='viridis', alpha=0.3)
    ax.clabel(contour, inline=True, fontsize=8)
    
    # Ograniczenia: x1 >= 1, x2 >= 1
    ax.axvline(x=1, color='red', linestyle='--', linewidth=1, alpha=0.5, label='x1 = 1')
    ax.axhline(y=1, color='red', linestyle='--', linewidth=1, alpha=0.5, label='x2 = 1')
    
    # Ograniczenie: sqrt(x1^2 + x2^2) <= a
    circle = Circle((0, 0), a_val, fill=False, edgecolor='red', linewidth=2, 
                    linestyle='-', label=f'r = {a_val}')
    ax.add_patch(circle)
    
    # Zaznaczenie rozwiązania dla zewnętrznej funkcji kary
    ax.plot(data['zewnetrzna']['x1'][idx], data['zewnetrzna']['x2'][idx], 
            'ro', markersize=10, label=f'Zewnętrzna: ({data["zewnetrzna"]["x1"][idx]:.2f}, {data["zewnetrzna"]["x2"][idx]:.2f})')
    
    # Zaznaczenie rozwiązania dla wewnętrznej funkcji kary
    ax.plot(data['wewnetrzna']['x1'][idx], data['wewnetrzna']['x2'][idx], 
            'bs', markersize=10, label=f'Wewnętrzna: ({data["wewnetrzna"]["x1"][idx]:.2f}, {data["wewnetrzna"]["x2"][idx]:.2f})')
    
    ax.set_xlabel('x1', fontsize=10)
    ax.set_ylabel('x2', fontsize=10)
    ax.set_title(f'a = {a_val}', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc='upper right')
    ax.set_xlim([0, 6])
    ax.set_ylim([0, 6])
    ax.set_aspect('equal')

plt.tight_layout()
plt.savefig('wyniki_optymalizacji_lab3.png', dpi=300, bbox_inches='tight')
plt.show()

# Dodatkowy wykres porównawczy
fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Wykres 1: Wartości funkcji celu vs parametr a
ax1.plot(data['a'], data['zewnetrzna']['y'], 'ro-', linewidth=2, markersize=8, label='Zewnętrzna funkcja kary')
ax1.plot(data['a'], data['wewnetrzna']['y'], 'bs-', linewidth=2, markersize=8, label='Wewnętrzna funkcja kary')
ax1.set_xlabel('Parametr a', fontsize=11)
ax1.set_ylabel('Wartość funkcji celu y*', fontsize=11)
ax1.set_title('Wartość funkcji celu w zależności od parametru a', fontsize=12, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=10)

# Wykres 2: Liczba wywołań funkcji celu vs parametr a
ax2.plot(data['a'], data['zewnetrzna']['calls'], 'ro-', linewidth=2, markersize=8, label='Zewnętrzna funkcja kary')
ax2.plot(data['a'], data['wewnetrzna']['calls'], 'bs-', linewidth=2, markersize=8, label='Wewnętrzna funkcja kary')
ax2.set_xlabel('Parametr a', fontsize=11)
ax2.set_ylabel('Liczba wywołań funkcji celu', fontsize=11)
ax2.set_title('Efektywność optymalizacji', fontsize=12, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.legend(fontsize=10)

plt.tight_layout()
plt.savefig('porownanie_wynikow_lab3.png', dpi=300, bbox_inches='tight')
plt.show()

# Wykres 3: Wszystkie rozwiązania na jednym wykresie
fig3, ax = plt.subplots(figsize=(10, 10))

# Wykres poziomic funkcji celu
contour = ax.contour(X1, X2, Z, levels=30, cmap='viridis', alpha=0.6)
contourf = ax.contourf(X1, X2, Z, levels=30, cmap='viridis', alpha=0.3)
ax.clabel(contour, inline=True, fontsize=8)

# Ograniczenia bazowe
ax.axvline(x=1, color='gray', linestyle='--', linewidth=1, alpha=0.5)
ax.axhline(y=1, color='gray', linestyle='--', linewidth=1, alpha=0.5)

# Wszystkie ograniczenia okręgiem
colors_circle = ['red', 'orange', 'purple']
for idx, a_val in enumerate(data['a']):
    circle = Circle((0, 0), a_val, fill=False, edgecolor=colors_circle[idx], 
                    linewidth=2, linestyle='-', label=f'a = {a_val}', alpha=0.7)
    ax.add_patch(circle)

# Wszystkie rozwiązania dla zewnętrznej funkcji kary
for idx, a_val in enumerate(data['a']):
    ax.plot(data['zewnetrzna']['x1'][idx], data['zewnetrzna']['x2'][idx], 
            'ro', markersize=10, alpha=0.7)
    ax.annotate(f'Z: a={a_val}', 
                (data['zewnetrzna']['x1'][idx], data['zewnetrzna']['x2'][idx]),
                xytext=(10, 10), textcoords='offset points', fontsize=9,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='red', alpha=0.3))

# Wszystkie rozwiązania dla wewnętrznej funkcji kary
for idx, a_val in enumerate(data['a']):
    ax.plot(data['wewnetrzna']['x1'][idx], data['wewnetrzna']['x2'][idx], 
            'bs', markersize=10, alpha=0.7)
    ax.annotate(f'W: a={a_val}', 
                (data['wewnetrzna']['x1'][idx], data['wewnetrzna']['x2'][idx]),
                xytext=(10, -20), textcoords='offset points', fontsize=9,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='blue', alpha=0.3))

ax.set_xlabel('x1', fontsize=12)
ax.set_ylabel('x2', fontsize=12)
ax.set_title('Wszystkie rozwiązania optymalizacji z funkcjami kary', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend(fontsize=10, loc='upper right')
ax.set_xlim([0, 6])
ax.set_ylim([0, 6])
ax.set_aspect('equal')

plt.tight_layout()
plt.savefig('wszystkie_rozwiazania_lab3.png', dpi=300, bbox_inches='tight')
plt.show()

print("Wykresy zostały wygenerowane i zapisane:")
print("1. wyniki_optymalizacji_lab3.png - trzy wykresy poziomic dla różnych wartości a")
print("2. porownanie_wynikow_lab3.png - porównanie wartości funkcji celu i liczby wywołań")
print("3. wszystkie_rozwiazania_lab3.png - wszystkie rozwiązania na jednym wykresie")
