"""
Numerical Simulation of the Two-Body Problem
Comparison of Euler Method and 4th Order Runge–Kutta Method

IB Mathematics: Analysis and Approaches – Internal Assessment
Candidate code: mjx609
"""

import numpy as np
import matplotlib.pyplot as plt
import csv

# =============================================================================
# USER-DEFINED SIMULATION PARAMETERS
# =============================================================================

T = 10.0
dt = float(input("Enter time step dt in years (e.g. 0.002): "))
steps = int(T / dt)
time = np.linspace(0, T, steps + 1)

# =============================================================================
# INITIAL CONDITIONS
# =============================================================================

x0, y0 = 1.0, 0.0
vx0, vy0 = 0.0, 2 * np.pi
R0 = np.sqrt(x0**2 + y0**2)
GM = 4 * np.pi**2

# =============================================================================
# NUMERICAL METHODS
# =============================================================================

def euler_step(x, y, vx, vy, dt):
    r = np.sqrt(x**2 + y**2)
    ax = -GM * x / r**3
    ay = -GM * y / r**3

    x_new = x + dt * vx
    y_new = y + dt * vy
    vx_new = vx + dt * ax
    vy_new = vy + dt * ay

    return x_new, y_new, vx_new, vy_new


def rk4_step(x, y, vx, vy, dt):
    def derivatives(state):
        x, y, vx, vy = state
        r = np.sqrt(x**2 + y**2)
        ax = -GM * x / r**3
        ay = -GM * y / r**3
        return np.array([vx, vy, ax, ay])

    state = np.array([x, y, vx, vy])

    k1 = derivatives(state)
    k2 = derivatives(state + dt * k1 / 2)
    k3 = derivatives(state + dt * k2 / 2)
    k4 = derivatives(state + dt * k3)

    return state + dt/6 * (k1 + 2*k2 + 2*k3 + k4)

# =============================================================================
# STORAGE
# =============================================================================

x_euler, y_euler = [x0], [y0]
x_rk4, y_rk4 = [x0], [y0]

radius_error_euler = []
radius_error_rk4 = []

absolute_radius_error_euler = []
absolute_radius_error_rk4 = []

theta_euler = [np.arctan2(y0, x0)]
theta_rk4 = [np.arctan2(y0, x0)]

# =============================================================================
# MAIN LOOP
# =============================================================================

x_e, y_e, vx_e, vy_e = x0, y0, vx0, vy0
x_r, y_r, vx_r, vy_r = x0, y0, vx0, vy0

for n in range(steps):
    x_e, y_e, vx_e, vy_e = euler_step(x_e, y_e, vx_e, vy_e, dt)
    x_euler.append(x_e)
    y_euler.append(y_e)
    theta_euler.append(np.arctan2(y_e, x_e))

    x_r, y_r, vx_r, vy_r = rk4_step(x_r, y_r, vx_r, vy_r, dt)
    x_rk4.append(x_r)
    y_rk4.append(y_r)
    theta_rk4.append(np.arctan2(y_r, x_r))

    r_e = np.sqrt(x_e**2 + y_e**2)
    r_r = np.sqrt(x_r**2 + y_r**2)

    radius_error_euler.append(abs(r_e - R0) / R0)
    radius_error_rk4.append(abs(r_r - R0) / R0)

    absolute_radius_error_euler.append(abs(r_e - 1.0))
    absolute_radius_error_rk4.append(abs(r_r - 1.0))

# =============================================================================
# FINAL RELATIVE RADIUS ERRORS
# =============================================================================

final_abs_error_euler = absolute_radius_error_euler[-1]
final_abs_error_rk4 = absolute_radius_error_rk4[-1]

print(f"Euler absolute radius error: {final_abs_error_euler:.3e} AU")
print(f"RK4 absolute radius error: {final_abs_error_rk4:.3e} AU")

final_relative_error_euler = radius_error_euler[-1]
final_relative_error_rk4 = radius_error_rk4[-1]

print(f"Euler relative radius error:   {final_relative_error_euler:.3e}")
print(f"RK4 relative radius error:     {final_relative_error_rk4:.3e}")
print(f"Difference in orders of magnitude: {np.log10(final_relative_error_euler / final_relative_error_rk4):.1f}")
print(f"RK4 is {final_relative_error_euler / final_relative_error_rk4:.1e} times more accurate")

# =============================================================================
# NUMBER OF REVOLUTIONS
# =============================================================================

theta_euler = np.arctan2(np.array(y_euler), np.array(x_euler))
theta_euler_unwrapped = np.unwrap(theta_euler)
revolutions_euler = theta_euler_unwrapped[-1] / (2 * np.pi)
revolutions_euler_rounded = round(revolutions_euler, 4)

print(f"Euler completed revolutions: {revolutions_euler_rounded}")

theta_rk4 = np.arctan2(np.array(y_rk4), np.array(x_rk4))
theta_rk4_unwrapped = np.unwrap(theta_rk4)
revolutions_rk4 = theta_rk4_unwrapped[-1] / (2 * np.pi)
revolutions_rk4_rounded = round(revolutions_rk4, 4)

print(f"RK4 completed revolutions: {revolutions_rk4_rounded}")

# =============================================================================
# CSV
# =============================================================================

with open("euler.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "time", "x", "y",
        "relative_radius_error", "absolute_radius_error"
    ])

    for i in range(0, steps, 10):
        writer.writerow([
            round(time[i], 4),
            round(x_euler[i], 6),
            round(y_euler[i], 6),
            f"{radius_error_euler[i-1]:.3e}" if i > 0 else "0",
            f"{absolute_radius_error_euler[i-1]:.3e}" if i > 0 else "0"
        ])

with open("rk4.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "time", "x", "y",
        "relative_radius_error", "absolute_radius_error"
    ])

    for i in range(0, steps, 10):
        writer.writerow([
            round(time[i], 4),
            round(x_rk4[i], 6),
            round(y_rk4[i], 6),
            f"{radius_error_rk4[i-1]:.3e}" if i > 0 else "0",
            f"{absolute_radius_error_rk4[i-1]:.3e}" if i > 0 else "0"
        ])

# =============================================================================
# VISUALISATION
# =============================================================================

theta = np.linspace(0, 2*np.pi, 400)

plt.figure(figsize=(8, 8))
plt.plot(np.cos(theta), np.sin(theta), 'k:', alpha=0.3, label="Exact Orbit")
plt.plot(x_rk4, y_rk4, label="RK4", linewidth=2)
plt.plot(x_euler, y_euler, "--", label="Euler")
plt.scatter(0, 0, color="orange", label="Sun")
plt.axis("equal")
plt.xlabel("x (AU)")
plt.ylabel("y (AU)")
plt.title("Orbital Trajectories")
plt.legend()
plt.show()

plt.figure(figsize=(10, 5))
plt.plot(time[1:], radius_error_euler, label="Euler Method")
plt.plot(time[1:], radius_error_rk4, label="Runge–Kutta 4")
plt.yscale("log")
plt.xlabel("Time (years)")
plt.ylabel("Relative Radius Error")
plt.title("Relative Radius Error Over Time")
plt.legend()
plt.grid(alpha=0.3)
plt.show()

plt.figure(figsize=(10, 5))
plt.plot(time[1:], absolute_radius_error_euler, label="Euler Method")
plt.plot(time[1:], absolute_radius_error_rk4, label="RK4 Method")
plt.xlabel("Time (years)")
plt.ylabel("Absolute Radius Error (AU)")
plt.title("Absolute Error of Orbital Radius")
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
