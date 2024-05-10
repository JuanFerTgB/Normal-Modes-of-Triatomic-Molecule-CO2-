import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle

# Physical constants
m = 16.0  # Mass of outer atoms (Oxygen)
M = 12.01  # Mass of central atom (Carbon)
k = 100.0  # Spring constant

# Frequencies and eigenvectors of normal modes
frequencies = [0, np.sqrt(k / m), np.sqrt((k * (2 * m + M)) / (m * M))]  # Frequencies of the normal modes
modes = [np.array([1, 1, 1]), np.array([1, 0, -1]), np.array([1, (-2 * m) / M, 1])]  # Eigenvectors of the normal modes

# Simulation configuration
dt = 0.05  # Time step
t_max = 10  # Maximum time
time = np.arange(0, t_max, dt)  # Time array

# Function to calculate atom positions
def calculate_positions(omega, mode, time):
    A = 0.17  # Amplitude
    positions = A * np.outer(np.cos(omega * time), mode)  # Calculate positions based on the normal mode
    if mode[1] == (-2 * m) / M:  # Specific adjustment for the second vibrational mode
        max_displacement = abs(positions[:, 1].max())  # Maximum displacement of central atom
        positions[:, 0] = np.clip(positions[:, 0], -max_displacement, max_displacement)  # Clip x-coordinate
        positions[:, 2] = np.clip(positions[:, 2], -max_displacement, max_displacement)  # Clip z-coordinate
    return positions + np.array([0, 1, 2])  # Adjust positions to avoid overlapping with the fixed central atom

# Function to draw a spring between two points
def draw_spring(ax, start, end, num_coils, width=0.05):
    start = np.array(start)
    end = np.array(end)
    length = np.linalg.norm(end - start)  # Length of the spring
    vec = (end - start) / length  # Unit vector along the spring
    perp_vec = np.array([-vec[1], vec[0]]) * width  # Perpendicular vector for width
    t = np.linspace(0, length, num_coils * 10)
    x = np.linspace(start[0], end[0], num_coils * 10)
    y = np.linspace(start[1], end[1], num_coils * 10)
    sine_wave = np.sin(t * np.pi * num_coils / length) * perp_vec[1]  # Sinusoidal wave for spring shape
    spring_x = x + sine_wave * vec[1]  # x-coordinates of spring
    spring_y = y - sine_wave * vec[0]  # y-coordinates of spring
    line, = ax.plot(spring_x, spring_y, 'k-')  # Plot the spring
    return line

# Create an animation for each mode
for idx, omega in enumerate(frequencies):
    fig, ax = plt.subplots()  # Create figure and axis
    ax.set_xlim(-1, 3)  # Set x-axis limits
    ax.set_ylim(-1, 1)  # Set y-axis limits
    ax.axis('off')  # Turn off axis
    ax.set_title(f"Normal Mode {idx + 1} with Frequency {omega:.2f} Hz", fontsize=14)  # Set title of the plot

    # Initialize atoms and springs
    colors = ['red' if i == 1 else 'gray' for i in range(3)]  # Colors for atoms
    atoms = [Circle((i, 0), 0.16 if i == 1 else 0.1, fc=colors[i], ec='black', linewidth=1.5, alpha=1, zorder=3) for i in range(3)]  # Create circles for atoms
    texts = ['' if i == 1 else '' for i in range(3)]  # Texts for atoms (empty for central atom)
    text_objects = []  # Text objects list

    for i, atom in enumerate(atoms):  # Loop over atoms
        ax.add_patch(atom)  # Add atom to plot
        txt = ax.text(atom.center[0], atom.center[1], texts[i], color='white', ha='center', va='center', fontsize=12, zorder=4)  # Add text to atom
        text_objects.append(txt)  # Append text object to list
    
    springs = []  # Springs list

    def init():
        global springs  # Global variable for springs
        springs = [draw_spring(ax, atoms[i].center, atoms[i+1].center, 10) for i in range(len(atoms)-1)]  # Draw springs between atoms
        return atoms + springs  # Return atoms and springs

    def animate(i):
        positions = calculate_positions(omega, modes[idx], time[i])  # Calculate atom positions for current time
        for j, atom in enumerate(atoms):  # Loop over atoms
            atom.center = (positions[0,j], 0)  # Update atom position
        for k in range(len(springs)):  # Loop over springs
            springs[k].remove()  # Remove previous spring
            springs[k] = draw_spring(ax, atoms[k].center, atoms[k+1].center, 10)  # Draw new spring
        return atoms + springs  # Return updated atoms and springs

    ani = animation.FuncAnimation(fig, animate, init_func=init, frames=len(time),  # Create animation
                                  interval=50, blit=False)  # Set animation parameters
    path = f"C:\\Users\\Enrique Riasgos\\OneDrive\\Escritorio\\SEPTIMO SEMESTRE\\Mecanica clasica\\Tarea 3\\Normal_Mode{idx + 1}.gif"  # Path to save animation
    try:
        ani.save(path, writer='pillow', fps=20)  # Save animation as GIF
    except ValueError:
        print("Using Pillow instead of ImageMagick.")  # Print message if Pillow is used as a writer
        ani.save(path, writer='pillow', fps=20)  # Save animation as GIF
    plt.close(fig)  # Close figure


