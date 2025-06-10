import sys
import csv
import numpy as np
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QWidget, QPushButton,
    QComboBox, QSlider, QLabel, QHBoxLayout, QCheckBox, QFileDialog
)
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT
from matplotlib.figure import Figure
import matplotlib.pyplot as plt


class PickGUI(QMainWindow):
    def __init__(self, cradar_obj):
        super().__init__()
        self.setWindowTitle("Radargram Picker")
        self.cradar = cradar_obj
        self.picking = False
        self.deleting = False
        self.picked_layers = {}  # Dictionary: layer_name -> list of segments
        self.current_layer_name = None
        self.segments = [[]]  # List of segments (each is list of (x, y))
        self.pick_lines = []  # Store pick line references

        self.init_ui()

    def init_layer_checkboxes(self):
        self.layer_checkboxes = {}
        for i in reversed(range(self.layer_checkbox_layout.count())):
            widget = self.layer_checkbox_layout.itemAt(i).widget()
            if widget is not None:
                widget.setParent(None)
        for layer_name in self.cradar.Layer.keys():
            checkbox = QCheckBox(layer_name)
            checkbox.setChecked(True)
            checkbox.stateChanged.connect(self.plot_layers)
            self.layer_checkboxes[layer_name] = checkbox
            self.layer_checkbox_layout.addWidget(checkbox)

    def init_ui(self):
        main_widget = QWidget()
        layout = QVBoxLayout(main_widget)

        self.canvas = FigureCanvas(Figure())
        self.canvas.setSizePolicy(self.canvas.sizePolicy().Expanding, self.canvas.sizePolicy().Expanding)
        self.ax = self.canvas.figure.subplots()
        self.im = None
        self.canvas.mpl_connect("button_press_event", self.on_click)

        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

        self.cmap_selector = QComboBox()
        self.cmap_selector.addItems(plt.colormaps())
        self.cmap_selector.setCurrentText("bone_r")
        self.cmap_selector.currentTextChanged.connect(self.update_colormap)

        data = self.cradar.Data
        if hasattr(data, "values"):
            data = data.values
        self.data_min = int(np.min(data))
        self.data_max = int(np.max(data))

        self.vmin_slider = QSlider(Qt.Horizontal)
        self.vmin_slider.setMinimum(self.data_min)
        self.vmin_slider.setMaximum(self.data_max)
        self.vmin_slider.setValue(self.data_min)
        self.vmin_slider.valueChanged.connect(self.update_colormap)

        self.vmax_slider = QSlider(Qt.Horizontal)
        self.vmax_slider.setMinimum(self.data_min)
        self.vmax_slider.setMaximum(self.data_max)
        self.vmax_slider.setValue(self.data_max)
        self.vmax_slider.valueChanged.connect(self.update_colormap)

        self.vmin_label = QLabel(f"vmin: {self.data_min}")
        self.vmax_label = QLabel(f"vmax: {self.data_max}")

        slider_layout = QHBoxLayout()
        slider_layout.addWidget(self.vmin_label)
        slider_layout.addWidget(self.vmin_slider)
        slider_layout.addWidget(self.vmax_label)
        slider_layout.addWidget(self.vmax_slider)

        self.pick_button = QCheckBox("Pick-Modus")
        self.pick_button.stateChanged.connect(self.toggle_picking)

        self.delete_button = QCheckBox("Delete-Modus")
        self.delete_button.stateChanged.connect(self.toggle_deleting)

        self.stop_segment_button = QPushButton("Stop Segment")
        self.stop_segment_button.clicked.connect(self.stop_segment)

        self.export_button = QPushButton("Export Picks")
        self.export_button.clicked.connect(self.export_picks)

        self.load_layer_button = QPushButton("Layer anzeigen")
        self.load_layer_button.clicked.connect(self.init_layer_checkboxes)

        self.layer_checkbox_container = QWidget()
        self.layer_checkbox_layout = QVBoxLayout()
        self.layer_checkbox_container.setLayout(self.layer_checkbox_layout)
        layout.addWidget(self.layer_checkbox_container)
        self.load_layer_button.clicked.connect(self.plot_layers)

        layout.addWidget(QLabel("Colormap"))
        layout.addWidget(self.cmap_selector)
        layout.addLayout(slider_layout)
        layout.addWidget(self.pick_button)
        layout.addWidget(self.delete_button)
        layout.addWidget(self.stop_segment_button)
        layout.addWidget(self.export_button)
        layout.addWidget(self.load_layer_button)

        self.setCentralWidget(main_widget)
        self.plot_radargram()

    def plot_radargram(self):
        # Speichere aktuellen View-Zustand
        try:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
        except:
            xlim, ylim = None, None
        data = self.cradar.Data
        if hasattr(data, "values"):
            data = data.values

        vmin = self.vmin_slider.value()
        vmax = self.vmax_slider.value()

        self.vmin_label.setText(f"vmin: {vmin}")
        self.vmax_label.setText(f"vmax: {vmax}")

        cmap = self.cmap_selector.currentText()

        if self.im is None:
            self.im = self.ax.imshow(data, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
            self.ax.set_title(f"{self.cradar.Frame}")
            self.ax.set_ylim(data.shape[0], 0)
            self.ax.set_xlim(0, data.shape[1])
        else:
            self.im.set_data(data)
            self.im.set_clim(vmin, vmax)
            self.im.set_cmap(cmap)

        # if xlim is not None and ylim is not None:
        #     self.ax.set_xlim(xlim)
        #     self.ax.set_ylim(ylim)

        self.canvas.draw()
        self.plot_picked_points()

    def update_colormap(self):
        self.plot_radargram()

    def toggle_picking(self):
        if self.pick_button.isChecked():
            options = list(self.cradar.Layer.keys()) + ["Neuer Horizont"]
            name_dialog = QComboBox()
            name_dialog.addItems(options)
            name_dialog.setEditable(True)
            name_dialog.setCurrentText("Neuer Horizont")

            name, _ = QFileDialog.getSaveFileName(self, "Horizontname wählen oder bestehenden auswählen", "", "")
            if name:
                self.current_layer_name = name
            else:
                self.current_layer_name = f"Horizont_{len(self.picked_layers)+1}"

            # Lade bestehende Daten aus self.picked_layers oder self.cradar.Layer
            if self.current_layer_name in self.picked_layers:
                self.segments = [list(seg) for seg in self.picked_layers[self.current_layer_name]]
            elif self.current_layer_name in self.cradar.Layer:
                layer = self.cradar.Layer[self.current_layer_name]
                x = layer["trace"]
                if "value_idx" in layer:
                    y = [idx if idx < len(self.cradar.Time) else np.nan for idx in layer["value_idx"]]
                else:
                    y = layer["value"]
                self.segments = [[(x[i], y[i]) for i in range(len(x))]]
            else:
                self.segments = [[]]

            self.picking = True
            self.delete_button.setChecked(False)
        else:
            self.picking = False
            self.current_layer_name = None

    def toggle_deleting(self):
        self.deleting = self.delete_button.isChecked()
        if self.deleting:
            self.pick_button.setChecked(False)

    def stop_segment(self):
        if self.segments and self.segments[-1]:
            self.picked_layers.setdefault(self.current_layer_name, []).append(self.segments[-1])
            self.segments.append([])

        if self.segments[-1]:
            self.segments.append([])

    def on_click(self, event):
        if not event.inaxes:
            return

        x = int(event.xdata)
        y = float(event.ydata)

        if self.deleting:
            tol = 10
            for segment in self.segments:
                for i, (px, py) in enumerate(segment):
                    if abs(px - x) < tol and abs(py - y) < tol:
                        segment.pop(i)
                        break
            self.plot_picked_points()

        elif self.picking:
            if not self.segments:
                self.segments = [[]]
            if not self.segments[-1]:
                for segment in self.segments:
                    if segment and (abs(segment[-1][0] - x) < 5 and abs(segment[-1][1] - y) < 5):
                        self.segments[-1] = segment
                        break
            self.segments[-1].append((x, y))
            self.plot_picked_points()

    def plot_picked_points(self):
        for line in self.pick_lines:
            line.remove()
        self.pick_lines = []

        for segment in self.segments:
            if segment:
                x_vals, y_vals = zip(*segment)
                line, = self.ax.plot(x_vals, y_vals, 'r-', marker='o')
                self.pick_lines.append(line)
        self.canvas.draw()

    def export_picks(self):
        from scipy.interpolate import interp1d
        for layer_name, segments in self.picked_layers.items():
            all_traces = []
            all_values = []
            for seg in segments:
                if len(seg) < 2:
                    continue
                x_vals, y_vals = zip(*seg)
                f = interp1d(x_vals, y_vals, bounds_error=False, fill_value="extrapolate")
                x_full = np.arange(min(x_vals), max(x_vals)+1)
                y_interp = f(x_full)
                all_traces.extend(x_full.tolist())
                all_values.extend(y_interp.tolist())
            self.cradar.add_layer_by_frame_trace(layer_name, all_traces, all_values, color=(227, 26, 28))

            # Exportiere jeden Layer als eigene Datei
            safe_name = layer_name.replace(" ", "_").replace("/", "_").replace("\\", "_")
            default_name = f"{str(self.cradar.Frame)}_{safe_name}.csv"
            filename, _ = QFileDialog.getSaveFileName(self, f"Speichern für {layer_name}", default_name, "CSV files (*.csv)")
            if filename:
                with open(filename, mode='w', newline='') as file:
                    writer = csv.writer(file, delimiter='	')
                    writer.writerow(["ProfileID", "TraceNumber", "TWT"])
                    for trace, val in zip(all_traces, all_values):
                        try:
                            twt = self.cradar.Time[int(val)] if int(val) < len(self.cradar.Time) else ''
                        except:
                            twt = ''
                        writer.writerow([self.cradar.Frame, trace, twt])
        print("==> Alle Picks exportiert")

    def plot_layers(self):
        # Entferne nur vorhandene Layer-Linien
        for line in self.ax.lines[:]:
            label = line.get_label()
            if label in self.layer_checkboxes:
                line.remove()

        # Zeichne nur die aktivierten Layer
        for layer_name, checkbox in self.layer_checkboxes.items():
            if checkbox.isChecked():
                layer = self.cradar.Layer[layer_name]
                if "value_idx" in layer:
                    x = layer["trace"]
                    y = [idx if idx < len(self.cradar.Time) else np.nan for idx in layer["value_idx"]]
                else:
                    x = layer["trace"]
                    y = layer["value"]

                # Segmentierung bei Trace-Lücken
                x_segments = []
                y_segments = []
                segment_x = [x[0]]
                segment_y = [y[0]]
                for i in range(1, len(x)):
                    if x[i] - x[i - 1] > 1:
                        x_segments.append(segment_x)
                        y_segments.append(segment_y)
                        segment_x = [x[i]]
                        segment_y = [y[i]]
                    else:
                        segment_x.append(x[i])
                        segment_y.append(y[i])
                x_segments.append(segment_x)
                y_segments.append(segment_y)

                for sx, sy in zip(x_segments, y_segments):
                    color = "red" if layer_name.lower() == "base" else "yellow" if layer_name.lower() == "surface" else None
                    self.ax.plot(sx, sy, label=layer_name, linestyle='-', marker='.', color=color, markersize=0.05)

        self.ax.legend()
        self.canvas.draw()