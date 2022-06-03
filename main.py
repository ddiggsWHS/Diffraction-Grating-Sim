from math import sin, cos, pi, atan, asin
import numpy as np
import matplotlib.pyplot as plt


def Intensity(wall_x, aperture_x):
    if (pi * a / (wavelength * D) * (wall_x - aperture_x)) ** 2 != 0:
        return I0 * sin(pi * a / (wavelength * D) * (wall_x - aperture_x)) ** 2 / (
                pi * a / (wavelength * D) * (wall_x - aperture_x)) ** 2
    else:
        return I0


def Phase(wall_x, aperture_x):
    return (2 * pi / wavelength * ((wall_x - aperture_x) ** 2 + D ** 2) ** 0.5) % (2 * pi)


def combine_sine_waves(phase1, A1, phase2, A2):
    phase_final = atan((A1 * sin(phase1) + A2 * sin(phase2)) / (A1 * cos(phase1) + A2 * cos(phase2)))
    amp_final = (A1 * sin(phase1) + A2 * sin(phase2)) / sin(phase_final)
    return amp_final, phase_final


def Wall_Plot(list_of_slits, min_wall_x, max_wall_x, wall_step_size):
    wall_points = np.arange(min_wall_x, max_wall_x, wall_step_size)
    diffraction_points = np.zeros(wall_points.shape)
    phases = np.empty((len(list_of_slits), wall_points.shape[0]))
    intensity_points = np.zeros(wall_points.shape)
    interference_points = np.zeros(wall_points.shape)
    phase_final = np.zeros(wall_points.shape)
    amp_final = np.zeros(wall_points.shape)
    for i, slit in enumerate(list_of_slits):
        for j, point in enumerate(wall_points):
            diffraction_points[j] += Intensity(point, slit)
            phases[i][j] = Phase(point, slit)
    phase_final[:] = phases[0][:]
    amp_final[:] = 1
    for i, slit in enumerate(list_of_slits):
        for j, point in enumerate(wall_points):
            if i != 0:
                temp_amp, temp_phase = combine_sine_waves(phase_final[j], amp_final[j], phases[i][j], 1)
                phase_final[j] = temp_phase
                amp_final[j] = temp_amp
    for j in range(len(wall_points)):
        interference_points[j] = amp_final[j] ** 2
    diffraction_points = diffraction_points / max(diffraction_points)
    interference_points = (interference_points - min(interference_points)) / max(interference_points)
    for i, point in enumerate(wall_points):
        intensity_points[i] = diffraction_points[i] * interference_points[i]
    intensity_points = abs(intensity_points)
    intensity_points = intensity_points / max(intensity_points)
    return wall_points, intensity_points, diffraction_points, interference_points


def main(angle, show_diffraction, show_interference, show_intensity):
    try:
        assert D / a > 10
        assert D / (a ** 2 / wavelength) > 10
    except AssertionError:
        print('Fraunhofer conditions NOT met')
        print('Both these values need to be above 10')
        print('D / a')
        print(D / a)
        print('D / (a^2 / wavelength)')
        print(D / (a ** 2 / wavelength))
        quit()
    else:
        print('Fraunhofer conditions met')
    slits = [-slit_separation / 2 + slit_separation * n for n in range(slit_count)]
    points, intensity, diffraction, interference = Wall_Plot(slits, -D, D, step_size)
    if angle:
        x = [atan(x / D) * 180 / pi for x in points]
        plt.xticks([-45, -30, -15, 0, 15, 30, 45])
    else:
        x = points
    fig, ax = plt.subplots(figsize=(10, 10))
    if show_diffraction:
        ax.plot(x, diffraction)
    if show_interference:
        ax.plot(x, interference)
    if show_intensity:
        ax.plot(x, intensity)
    plt.savefig('{}_slit.svg'.format(slit_count))
    print('done')


wavelength = 660 * 10 ** -9
D = 10
I0 = 1
a = 2 * wavelength  # slit thickness
step_size = 2 * D / 10 ** 5
slit_count = 100
slit_separation = 6 * wavelength

main(False, False, False, True)
