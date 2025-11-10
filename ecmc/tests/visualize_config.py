#!/usr/bin/env python3
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import sys

def read_configuration(filename):
    """Read configuration file
    
    Format:
        Line 1: Lx Ly
        Following lines: x y radius
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # First line: box size
    Lx, Ly = map(float, lines[0].split())
    
    # Read particle data
    particles = []
    for line in lines[1:]:
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue
        x, y, r = map(float, line.split())
        particles.append((x, y, r))
    
    return Lx, Ly, particles

def plot_configuration(filename, title="Particle Configuration", save_as=None):
    """Plot configuration"""
    Lx, Ly, particles = read_configuration(filename)
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    
    # Set axis range (centered box: [-L/2, L/2])
    ax.set_xlim(-Lx/2, Lx/2)
    ax.set_ylim(-Ly/2, Ly/2)
    ax.set_aspect('equal')
    
    # Draw particles
    for x, y, r in particles:
        circle = patches.Circle((x, y), r, 
                                edgecolor='black', 
                                facecolor='lightblue', 
                                linewidth=1.5,
                                alpha=0.7)
        ax.add_patch(circle)
    
    # Draw box boundary (centered)
    box = patches.Rectangle((-Lx/2, -Ly/2), Lx, Ly, 
                            linewidth=2, 
                            edgecolor='red', 
                            facecolor='none')
    ax.add_patch(box)
    
    ax.set_xlabel('X', fontsize=14)
    ax.set_ylabel('Y', fontsize=14)
    ax.set_title(f'{title}\nN={len(particles)}, Box=[{Lx:.1f}, {Ly:.1f}] (centered at origin)', 
                 fontsize=16)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    
    plt.tight_layout()
    
    if save_as:
        plt.savefig(save_as, dpi=150, bbox_inches='tight')
        print(f"Figure saved to {save_as}")
    
    plt.show()

def compare_configurations(file1, file2, title1="Initial", title2="Final", save_as=None):
    """Compare two configurations"""
    Lx1, Ly1, particles1 = read_configuration(file1)
    Lx2, Ly2, particles2 = read_configuration(file2)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    for ax, Lx, Ly, particles, title in [
        (ax1, Lx1, Ly1, particles1, title1),
        (ax2, Lx2, Ly2, particles2, title2)
    ]:
        # Set axis range (centered box: [-L/2, L/2])
        ax.set_xlim(-Lx/2, Lx/2)
        ax.set_ylim(-Ly/2, Ly/2)
        ax.set_aspect('equal')
        
        # Draw particles
        for x, y, r in particles:
            circle = patches.Circle((x, y), r, 
                                    edgecolor='black', 
                                    facecolor='lightblue', 
                                    linewidth=1.5,
                                    alpha=0.7)
            ax.add_patch(circle)
        
        # Draw box boundary (centered)
        box = patches.Rectangle((-Lx/2, -Ly/2), Lx, Ly, 
                                linewidth=2, 
                                edgecolor='red', 
                                facecolor='none')
        ax.add_patch(box)
        
        ax.set_xlabel('X', fontsize=12)
        ax.set_ylabel('Y', fontsize=12)
        ax.set_title(f'{title}\nN={len(particles)}', fontsize=14)
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
        ax.axvline(x=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    
    plt.tight_layout()
    
    if save_as:
        plt.savefig(save_as, dpi=150, bbox_inches='tight')
        print(f"Comparison saved to {save_as}")
    
    plt.show()

if __name__ == "__main__":
    import os
    
    # Try to find configuration files in multiple locations
    search_dirs = [
        ".",                    # Current directory
        "build",                # build subdirectory
        "../build",             # Parent directory's build
        os.path.dirname(__file__) + "/build",  # Script directory's build
    ]
    
    initial_file = None
    final_file = None
    
    for dir_path in search_dirs:
        init_path = os.path.join(dir_path, "config_initial.dat")
        final_path = os.path.join(dir_path, "config_final.dat")
        
        if os.path.exists(init_path):
            initial_file = init_path
        if os.path.exists(final_path):
            final_file = final_path
        
        if initial_file and final_file:
            break
    
    # Check if files exist
    if initial_file and final_file:
        print(f"Found both configurations in: {os.path.dirname(initial_file)}")
        print("Plotting comparison...")
        compare_configurations(initial_file, final_file,
                              title1="Initial Configuration",
                              title2="Final Configuration (after ECMC)",
                              save_as="comparison.png")
    elif initial_file:
        print(f"Found initial configuration in: {os.path.dirname(initial_file)}")
        plot_configuration(initial_file, 
                          title="Initial Configuration",
                          save_as="initial.png")
    elif final_file:
        print(f"Found final configuration in: {os.path.dirname(final_file)}")
        plot_configuration(final_file, 
                          title="Final Configuration",
                          save_as="final.png")
    else:
        print("No configuration files found!")
        print(f"Searched in: {search_dirs}")
        print("Please run the test_high_density executable first.")
        sys.exit(1)
