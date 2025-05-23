import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import numpy as np

# Define work packages with start and end months
work_packages = {
    "WP1": (1, 6),
    "WP2": (3, 46),
    "WP3": (10, 30),
    "WP4": (12, 48)
}

# Define milestones: (month, label, associated WP)
milestones = [
    (6, "MS1", "WP1"),
    (18, "MS2", "WP2"),
    (24, "MS3", "WP3"),
    (43, "MS4", "WP4")
]

# Setup
fig, ax = plt.subplots(figsize=(10, 3.2))
cmap = get_cmap("Dark2")
wp_names = list(work_packages.keys())
y_positions = np.arange(len(wp_names))[::-1]  # Top-down ordering
colors = [cmap(i) for i in range(len(wp_names))]

# Plot WP lines
for i, wp in enumerate(wp_names):
    start, end = work_packages[wp]
    ax.hlines(y=y_positions[i], xmin=start, xmax=end, color=colors[i], linewidth=10)

# Plot milestones with larger dots and bold labels
for month, label, wp in milestones:
    y = y_positions[wp_names.index(wp)]
    ax.plot(month, y, 'o', color=colors[wp_names.index(wp)], markersize=24)
    ax.text(month + 1, y + 0.2, label, fontsize=16, fontweight='bold', verticalalignment='center')

# Y-axis formatting (bold, larger font)
ax.set_yticks(y_positions)
ax.set_yticklabels(wp_names, fontsize=14, fontweight='bold')

# X-axis and title
ax.set_xticks(np.arange(1, 49, 3))
ax.set_xlim(0, 49)
ax.set_ylim(-0.75, len(wp_names) - 0.25)
ax.set_xlabel("Month", fontsize=12)
ax.set_title("Tensor Omics – Gantt Chart", fontsize=14)
ax.grid(True, axis='x', linestyle='--', linewidth=0.6)

plt.tight_layout()

# Save the figure
plt.savefig("tensor_omics_gantt-chart.pdf", format="pdf", dpi=300)
plt.savefig("tensor_omics_gantt-chart.png", format="png", dpi=300)

# Optional:
plt.show()
