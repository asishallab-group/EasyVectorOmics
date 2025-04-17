from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline
from sys import argv
np.random.seed(2025)
plt.rcParams.update({"figure.dpi": 200, "font.size": 15})
RELEVANT_ARROW_SIZE = 7
DOI_COLOR = "#003CE0"


def main():
    plot_1d()
    plt.savefig(argv[0].rsplit(".", 1)[0] + ".png", bbox_inches='tight')


def draw_smooth_line(ax, points, k=2, periodic=False, label="", label_move_rel=(0, 0), color="blue"):
    # Separate into x and y arrays
    x, y = zip(*points)

    # Parameterize the points
    t = np.linspace(0, 1, len(x))  # Parameterize as a sequence
    t_smooth = np.linspace(0, 1, 300)  # Dense parameterization for smoothing

    # Create smooth interpolation
    x_spline = make_interp_spline(t, x, k=k, bc_type="periodic" if periodic else None)  # Cubic interpolation
    y_spline = make_interp_spline(t, y, k=k, bc_type="periodic" if periodic else None)
    x_smooth = x_spline(t_smooth)
    y_smooth = y_spline(t_smooth)

    # Plot the smooth line
    ax.plot(x_smooth, y_smooth, linestyle='-', color=color)

    mid_index = len(x_smooth) // 2
    ax.text(x_smooth[mid_index] + label_move_rel[0], y_smooth[mid_index] + label_move_rel[1], label, color=color, ha="center")
    ax.set_aspect('equal')  # Maintain aspect ratio for shapes like circles


# refers to plot: https://docs.google.com/presentation/d/1T9e2uZamu__GOxcbz7fM81QUgAlvKp5JjUTwnGWbz3w/edit?slide=id.g349a373bb74_0_11#slide=id.g349a373bb74_0_11
def plot_1d():
    fig, ax = plt.subplots()
    ax.set_xlim(0, 125)
    ax.set_ylim(0, 110)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', length=0)
    xmin, xmax = ax.get_xlim()
    ax.set_xticks(list(range(int(xmin), int(xmax), 7)))
    ymin, ymax = ax.get_ylim()
    ax.set_yticks(list(range(int(ymin), int(ymax), 7)))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlabel(r"$tissue\ X$")
    ax.set_ylabel(r"$tissue\ Y$")
    ax.grid(True, linewidth=0.5, color='gray', alpha=0.2)
    add_compass(fig)

    # Add arrows at the ends of the axes
    ax.plot(1, 0, ">k", transform=ax.get_yaxis_transform(), clip_on=False)  # Arrow for x-axis
    ax.plot(0, 1, "^k", transform=ax.get_xaxis_transform(), clip_on=False)  # Arrow for y-axis

    def arrow(start_x, start_y, end_x, direction=1, length=12):
        end_y = start_y + direction * (length**2 - (abs(start_x) - abs(end_x))**2)**.5
        return ((start_x, start_y), (end_x, end_y))

    def standalone_arrows():
        def random_arrow(start_x, start_y, end_x, direction=1, length=None):
            if length is None:
                min_length = max(5, abs(start_x - end_x))
                length = np.random.uniform(min_length, RELEVANT_ARROW_SIZE)

            return arrow(start_x, start_y, end_x, direction, length=length)

        high_impact = (
            arrow(18, 56.5, 7, -1, length=14),
            arrow(10, 34.5, 17, length=13),
            arrow(25, 49.5, 22, -1, length=11),
        )
        return (
            *high_impact,
            random_arrow(7, 65, 9),
            random_arrow(112, 45, 113),
            random_arrow(90.5, 97.5, 97.5),
            random_arrow(21.5, 20, 15.5, -1),
            random_arrow(35.5, 105, 32, -1, length=7),
            random_arrow(35.5, 30, 35.5, -1),
            random_arrow(70, 40, 71, length=7),
            random_arrow(97.5, 5, 104.5),
            random_arrow(84, 12, 86, -1, length=6),
            random_arrow(49, 74, 48, -1),
            random_arrow(87, 84.5, 91, -1),
            random_arrow(84.5, 61, 80, -1),
            random_arrow(31, 6.5, 35),
            random_arrow(108, 21.5, 111, -1),
            random_arrow(65, 70.5, 68, -1),
            random_arrow(56, 40.5, 59, -1),
            random_arrow(117, 70.5, 120, -1),
            random_arrow(115, 28, 110, -1),
        )

    # first from left
    def left_blue(color="blue"):
        # Define the arrows, starting from top
        arrows = (
            arrow(21.5, 100, 17, length=10),
            arrow(14.5, 95, 7, length=11),
            arrow(21.5, 90, 16, length=10),
            arrow(28.5, 83, 23, length=11),
            arrow(28.5, 72, 19, length=13),
            arrow(42.5, 63, 35, length=13),
            arrow(35.5, 56, 31, length=10),
            arrow(44, 48.5, 38, length=14),
        )

        route_left_line = [
            (arrows[1][1][0] - 5, arrows[1][1][1] + 1),
            (arrows[4][1][0] - 5, arrows[4][1][1]),
            (arrows[6][1][0] - 5, arrows[6][1][1]),
            (arrows[7][0][0] - 9, arrows[7][0][1]),
        ]
        draw_smooth_line(ax, route_left_line, color=color, k=2)

        route_right_line = [
            (arrows[0][1][0] + 4, arrows[0][1][1] + 1),
            (arrows[3][0][0] + 5, arrows[3][0][1]),
            (arrows[5][0][0] + 3, arrows[5][0][1] + 3),
            (arrows[7][0][0] + 10, arrows[7][0][1]),
        ]
        draw_smooth_line(ax, route_right_line, color=color, k=2)

        return arrows

    # second from left, the one below the circle with two arrows
    def right_blue(color="blue"):
        arrows = (
            arrow(56, 16, 52),
            arrow(56, 10, 47),
            arrow(63, 2, 59, length=12)
        )

        route_left_line = [
            (arrows[0][1][0] - 10, arrows[0][1][1] - 2),
            (arrows[1][1][0] - 3, arrows[1][1][1] + 1),
            (arrows[1][0][0] - 7, arrows[1][0][1]),
            (arrows[1][0][0] - 3, arrows[1][0][1] - 4),
            (arrows[2][0][0] - 7, arrows[2][0][1]),
        ]
        draw_smooth_line(ax, route_left_line, color=color, label="DOI\nsubfields", label_move_rel=(-2, 22.5))

        route_right_line = [
            (arrows[0][1][0] + 5, arrows[0][1][1]),
            (arrows[0][0][0] + 5, arrows[0][0][1] + 1),
            (arrows[2][1][0] + 5, arrows[2][1][1] - 1),
            (arrows[2][0][0] + 5, arrows[2][0][1] + 3),
        ]
        draw_smooth_line(ax, route_right_line, color=color)

        return arrows

    # the circle
    def chaotic(color="blue"):
        arrows = (
            # upper
            arrow(63, 95, 64, length=10),

            # left from upper
            arrow(57, 98.5, 59, -1, length=9),

            # points to bottom-left
            arrow(63.5, 91, 55, -1, length=11),

            # bottom, points up
            arrow(59, 76.5, 62, length=8),

            # middle one, points down (not bottom right)
            arrow(65, 91.5, 66, -1, length=10),

            # middle one, points bottom-right
            arrow(66, 98.5, 71, -1, length=10),

            # points right, bottom one
            arrow(69.5, 86, 78, length=10),

            # points right, top one
            arrow(69.5, 98.5, 79.4, length=10),
        )
        # Define points for plotting
        route = [
            (start := (arrows[0][1][0] + 4, arrows[0][1][1] + 4)),
            (arrows[1][0][0] - 3, arrows[1][0][1]),
            (arrows[2][1][0] - 3, arrows[2][1][1] - 3),
            (arrows[3][0][0] + 3, arrows[3][0][1] - 1),
            (arrows[-2][1][0] + 3, arrows[-2][1][1] - 1),
            (arrows[-1][1][0] + 1, arrows[-1][1][1] + 4),
            start,
        ]

        # Combine points and make it periodic (ensure looping by appending the start point)
        draw_smooth_line(ax, route, periodic=True, color=color, label="Photosynthesis", label_move_rel=(0, 35))

        return arrows

    # the last one, first from right
    def merged(color="blue"):
        # Define the arrows, starting from top right
        merged = (
            arrow(105, 83.5, 116, length=14),
            arrow(101, 76, 102, length=10),
            arrow(94, 62.5, 100, length=13),
            arrow(90.5, 44, 100, length=23),
            arrow(97.5, 51, 105, length=12),
        )
        left_arm = (
            arrow(80, 27.5, 86, length=12),
            arrow(77, 15, 82, length=10),
        )
        right_arm = (
            arrow(95, 27, 95.5, length=10),
            arrow(98, 10, 94, length=14),
        )

        route_left_line = [
            (merged[0][1][0], merged[0][1][1] + 5),
            (merged[1][1][0] - 3, merged[1][1][1] + 1),
            (merged[2][1][0] - 6, merged[2][1][1]),
            (merged[2][0][0] - 4, merged[2][0][1]),

            (left_arm[0][1][0] - 5, left_arm[0][1][1]),
            (left_arm[1][0][0] - 6, left_arm[1][0][1]),
        ]
        draw_smooth_line(ax, route_left_line, color=color, k=2)

        route_right_line = [
            (merged[0][1][0] + 3, merged[0][1][1] - 5),
            (merged[1][0][0] + 12, merged[1][0][1]),
            (merged[4][1][0] + 3, merged[4][1][1] + 1),

            (right_arm[0][1][0] + 5, right_arm[0][1][1]),
            (right_arm[1][0][0] + 6, right_arm[1][0][1]),
        ]
        draw_smooth_line(ax, route_right_line, color=color, k=2, label="laminar\nflow", label_move_rel=(7, -29))

        route_bottom_line = [
            (left_arm[1][0][0] + 6, left_arm[1][0][1] - 3),
            (right_arm[1][1][0] - 6, right_arm[1][1][1]),
            (right_arm[1][0][0] - 6, right_arm[1][0][1]),
        ]
        draw_smooth_line(ax, route_bottom_line, color=color, k=2)
        return *merged, *left_arm, *right_arm

    arrows = (
        *standalone_arrows(),
        *left_blue(DOI_COLOR),
        *right_blue(DOI_COLOR),
        *chaotic("#7DB00F"),
        *merged("#D9591B")
    )
    for start, end in arrows:
        arrow_len = ((start[0] - end[0]) ** 2 + (start[1] - end[1]) ** 2) ** .5
        color = "black" if arrow_len > RELEVANT_ARROW_SIZE else "gray"
        ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle='->', color=color))


def add_compass(fig):
    # Add the inset axis for the compass, positioned to the left of the Y-axis
    compass_ax = fig.add_axes([0.704, 0.82, 0.1, 0.1], polar=True)  # Adjust position
    compass_ax.set_xticks([])
    compass_ax.set_yticks([])
    compass_ax.set_frame_on(False)

    # Add direction of interest (DOI)
    compass_ax.annotate('DOI', xy=(0, 0), xytext=(3*np.pi/4, 2),
                        arrowprops=dict(arrowstyle='<-', color=DOI_COLOR),
                        ha='center', color=DOI_COLOR)

    # Add compass arrows using annotate
    # Draw East arrow
    compass_ax.annotate('', xy=(0, 0.75), xytext=(0, 0),
                        arrowprops=dict(facecolor='black', shrink=0, width=1, headwidth=9, headlength=17.5),
                        ha='center', va='center')

    # Draw North arrow
    compass_ax.annotate('', xy=(np.pi/2, 0.75), xytext=(np.pi/2, 0),
                        arrowprops=dict(facecolor='black', shrink=0, width=10, headwidth=9, headlength=17.5),
                        ha='center', va='center')

    # Draw West arrow
    compass_ax.annotate('', xy=(np.pi, 0.75), xytext=(np.pi, 0),
                        arrowprops=dict(facecolor='black', shrink=0, width=1, headwidth=9, headlength=17.5),
                        ha='center', va='center')

    # Draw South arrow
    compass_ax.annotate('', xy=(3*np.pi/2, 0.75), xytext=(3*np.pi/2, 0),
                        arrowprops=dict(facecolor='black', shrink=0, width=1, headwidth=9, headlength=17.5),
                        ha='center', va='center')

    # Label directions
    directions = ['E', 'N', 'W', 'S']
    angles = [0, np.pi/2, np.pi, 3*np.pi/2]
    for angle, direction in zip(angles, directions):
        compass_ax.text(angle, 1.2, direction, ha='center', va='center', fontsize=13, fontweight='bold')


if __name__ == '__main__':
    main()
