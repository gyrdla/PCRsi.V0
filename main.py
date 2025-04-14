import sys
import json
import random
import numpy as np
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QTextEdit, QPushButton,
    QLabel, QFileDialog, QGroupBox, QRadioButton, QCheckBox, QLineEdit, QComboBox, QTabWidget,
    QMessageBox, QSizePolicy
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont, QPalette, QColor
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class PCRSimulator(QMainWindow):
    def __init__(self):
        super().__init__()
        self.dna_sequence = ""
        self.primers = {'forward': "", 'reverse': ""}
        self.probes = []
        self.db = self.load_database()
        self.init_ui()
        self.apply_styles()

    def init_ui(self):
        self.setWindowTitle("Clinical PCR Simulator")
        self.setFixedSize(800, 600)
        
        # Central widget with tabs
        self.tabs = QTabWidget()
        self.tabs.setTabPosition(QTabWidget.North)
        
        # Create tabs
        self.input_tab = self.create_input_tab()
        self.settings_tab = self.create_settings_tab()
        self.results_tab = self.create_results_tab()
        
        self.tabs.addTab(self.input_tab, "🧬 Input")
        self.tabs.addTab(self.settings_tab, "🔧 Settings")
        self.tabs.addTab(self.results_tab, "📊 Results")
        
        self.setCentralWidget(self.tabs)
        
        # Add control buttons
        self.create_control_buttons()

    def create_control_buttons(self):
        control_widget = QWidget()
        control_layout = QHBoxLayout()
        
        self.clear_btn = QPushButton("Clear")
        self.clear_btn.clicked.connect(self.clear_fields)
        
        self.sample_btn = QPushButton("Load Example")
        self.sample_btn.clicked.connect(self.load_example)
        
        self.run_btn = QPushButton("Run Simulation")
        self.run_btn.clicked.connect(self.run_simulation)
        
        control_layout.addWidget(self.clear_btn)
        control_layout.addWidget(self.sample_btn)
        control_layout.addWidget(self.run_btn)
        
        control_widget.setLayout(control_layout)
        
        # Add to the bottom of the window
        self.statusBar().addPermanentWidget(control_widget)

    def create_input_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        
        # DNA Sequence Input
        dna_group = QGroupBox("DNA Sequence")
        dna_layout = QVBoxLayout()
        self.dna_input = QTextEdit()
        self.dna_input.setPlaceholderText("Paste DNA sequence here or load from file...")
        self.load_fasta_btn = QPushButton("Load FASTA File")
        self.load_fasta_btn.clicked.connect(self.load_fasta)
        dna_layout.addWidget(self.dna_input)
        dna_layout.addWidget(self.load_fasta_btn)
        dna_group.setLayout(dna_layout)
        
        # Primer Input
        primer_group = QGroupBox("Primers")
        primer_layout = QVBoxLayout()
        
        forward_layout = QHBoxLayout()
        forward_layout.addWidget(QLabel("Forward:"))
        self.forward_primer = QLineEdit()
        self.forward_primer.setPlaceholderText("5'-3' sequence")
        forward_layout.addWidget(self.forward_primer)
        
        reverse_layout = QHBoxLayout()
        reverse_layout.addWidget(QLabel("Reverse:"))
        self.reverse_primer = QLineEdit()
        self.reverse_primer.setPlaceholderText("5'-3' sequence")
        reverse_layout.addWidget(self.reverse_primer)
        
        primer_layout.addLayout(forward_layout)
        primer_layout.addLayout(reverse_layout)
        primer_group.setLayout(primer_layout)
        
        # Probe Input
        probe_group = QGroupBox("Probe")
        probe_layout = QHBoxLayout()
        
        probe_layout.addWidget(QLabel("Sequence:"))
        self.probe_sequence = QLineEdit()
        self.probe_sequence.setPlaceholderText("5'-3' sequence")
        probe_layout.addWidget(self.probe_sequence)
        
        probe_layout.addWidget(QLabel("Dye:"))
        self.dye_combo = QComboBox()
        self.dye_combo.addItems(["FAM", "HEX", "VIC", "ROX", "CY5"])
        probe_layout.addWidget(self.dye_combo)
        
        self.multiplex_check = QCheckBox("Multiplex Mode")
        probe_layout.addWidget(self.multiplex_check)
        
        probe_group.setLayout(probe_layout)
        
        # Assemble tab
        layout.addWidget(dna_group)
        layout.addWidget(primer_group)
        layout.addWidget(probe_group)
        layout.addStretch()
        
        tab.setLayout(layout)
        return tab

    def create_settings_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        
        # PCR Settings
        pcr_group = QGroupBox("PCR Parameters")
        pcr_layout = QVBoxLayout()
        
        cycles_layout = QHBoxLayout()
        cycles_layout.addWidget(QLabel("Cycles:"))
        self.cycles_input = QLineEdit("40")
        cycles_layout.addWidget(self.cycles_input)
        pcr_layout.addLayout(cycles_layout)
        
        efficiency_layout = QHBoxLayout()
        efficiency_layout.addWidget(QLabel("Efficiency (0-1):"))
        self.efficiency_input = QLineEdit("0.95")
        efficiency_layout.addWidget(self.efficiency_input)
        pcr_layout.addLayout(efficiency_layout)
        
        threshold_layout = QHBoxLayout()
        threshold_layout.addWidget(QLabel("Threshold:"))
        self.threshold_input = QLineEdit("0.1")
        threshold_layout.addWidget(self.threshold_input)
        pcr_layout.addLayout(threshold_layout)
        
        pcr_group.setLayout(pcr_layout)
        
        # Thermal Profile
        thermal_group = QGroupBox("Thermal Profile")
        thermal_layout = QVBoxLayout()
        
        denature_layout = QHBoxLayout()
        denature_layout.addWidget(QLabel("Denaturation (°C):"))
        self.denature_temp = QLineEdit("95")
        denature_layout.addWidget(self.denature_temp)
        thermal_layout.addLayout(denature_layout)
        
        anneal_layout = QHBoxLayout()
        anneal_layout.addWidget(QLabel("Annealing (°C):"))
        self.anneal_temp = QLineEdit("60")
        anneal_layout.addWidget(self.anneal_temp)
        thermal_layout.addLayout(anneal_layout)
        
        extend_layout = QHBoxLayout()
        extend_layout.addWidget(QLabel("Extension (°C):"))
        self.extend_temp = QLineEdit("72")
        extend_layout.addWidget(self.extend_temp)
        thermal_layout.addLayout(extend_layout)
        
        thermal_group.setLayout(thermal_layout)
        
        # Assemble tab
        layout.addWidget(pcr_group)
        layout.addWidget(thermal_group)
        layout.addStretch()
        
        tab.setLayout(layout)
        return tab

    def create_results_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        
        # Results display
        results_group = QGroupBox("Simulation Results")
        results_layout = QVBoxLayout()
        
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111)
        
        results_layout.addWidget(self.canvas)
        
        # Results text
        self.results_text = QTextEdit()
        self.results_text.setReadOnly(True)
        results_layout.addWidget(self.results_text)
        
        results_group.setLayout(results_layout)
        
        layout.addWidget(results_group)
        tab.setLayout(layout)
        return tab

    def apply_styles(self):
        # Set font
        font = QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(10)
        self.setFont(font)
        
        # Set palette for light theme
        palette = QPalette()
        palette.setColor(QPalette.Window, QColor(240, 240, 240))
        palette.setColor(QPalette.WindowText, Qt.black)
        palette.setColor(QPalette.Base, QColor(255, 255, 255))
        palette.setColor(QPalette.AlternateBase, QColor(230, 230, 230))
        palette.setColor(QPalette.ToolTipBase, Qt.white)
        palette.setColor(QPalette.ToolTipText, Qt.black)
        palette.setColor(QPalette.Text, Qt.black)
        palette.setColor(QPalette.Button, QColor(240, 240, 240))
        palette.setColor(QPalette.ButtonText, Qt.black)
        palette.setColor(QPalette.BrightText, Qt.red)
        palette.setColor(QPalette.Highlight, QColor(0, 120, 215))
        palette.setColor(QPalette.HighlightedText, Qt.white)
        self.setPalette(palette)
        
        # Button styles
        button_style = """
        QPushButton {
            background-color: #f0f0f0;
            border: 1px solid #ccc;
            border-radius: 4px;
            padding: 5px 10px;
        }
        QPushButton:hover {
            background-color: #e0e0e0;
        }
        QPushButton:pressed {
            background-color: #d0d0d0;
        }
        """
        self.setStyleSheet(button_style)

    def load_database(self):
        try:
            with open('primer_probe_db.json') as file:
                return json.load(file)
        except (FileNotFoundError, json.JSONDecodeError):
            return {
                'targets': {
                    'SARS-CoV-2': {
                        'sequence': 'ATGTCTGATAATGGACCCCAAAATCAGCGAAATGCACCCCGCATTACGTTTGGTGGACCCTCAGATTCAACTGGCAGTAACCAGAATGGAGAACGCAGTGGGGCGCGATCAAAACAACGTCGGCCCCAAGGTTTACCCAATAATACTGCGTCTTGGTTCACCGCTCTCACTCAACATGGCAAGGAAGACCTTAAATTCCCTCGAGGACAAGGCGTTCCAATTAACACCAATAGCAGTCCAGATGACCAAATTGGCTACTACCGAAGAGCTACCAGACGAATTCGTGGTGGTGACGGTAAAATGAAAGATCTCAGTCCAAGATGGTATTTCTACTACCTAGGAACTGGGCCAGAAGCTGGACTTCCCTATGGTGCTAACAAAGACGGCATCATATGGGTTGCAACTGAGGGAGCCTTGAATACACCAAAAGATCACATTGGCACCCGCAATCCTGCTAACAATGCTGCAATCGTGCTACAACTTCCTCAAGGAACAACATTGCCAAAAGGCTTCTACGCAGAAGGGAGCAGAGGCGGCAGTCAAGCCTCTTCTCGTTCCTCATCACGTAGTCGCAACAGTTCAAGAAATTCAACTCCAGGCAGCAGTAGGGGAACTTCTCCTGCTAGAATGGCTGGCAATGGCGGTGATGCTGCTCTTGCTTTGCTGCTGCTTGACAGATTGAACCAGCTTGAGAGCAAAATGTCTGGTAAAGGCCAACAACAACAAGGCCAAACTGTCACTAAGAAATCTGCTGCTGAGGCTTCTAAGAAGCCTCGGCAAAAACGTACTGCCACTAAAGCATACAATGTAACACAAGCTTTCGGCAGACGTGGTCCAGAACAAACCCAAGGAAATTTTGGGGACCAGGAACTAATCAGACAAGGAACTGATTACAAACATTGGCCGCAAATTGCACAATTTGCCCCCAGCGCTTCAGCGTTCTTCGGAATGTCGCGCATTGGCATGGAAGTCACACCTTCGGGAACGTGGTTGACCTACACAGGTGCCATCAAATTGGATGACAAAGATCCAAATTTCAAAGATCAAGTCATTTTGCTGAATAAGCATATTGACGCATACAAAACATTCCCACCAACAGAGCCTAAAAAGGACAAAAAGAAGAAGGCTGATGAAACTCAAGCCTTACCGCAGAGACAGAAGAAACAGCAAACTGTGACTCTTCTTCCTGCTGCAGATTTGGATGATTTCTCCAAACAATTGCAACAATCCATGAGCAGTGCTGACTCAACTCAGGCCTAA',
                        'primers': {'forward': 'GACCCCAAAATCAGCGAAAT', 'reverse': 'TCTAGTTGCTGCTGCTACTG'},
                        'probe': {'sequence': 'ACCCCGCATTACGTTTGGTGG', 'dye': 'FAM', 'quencher': 'BHQ1'}
                    },
                    'HPV-16': {
                        'sequence': 'ATGCAACAGTTATGCACAGAGCTGCAAACAACTGCAATGTGTAGAGTTTCTGCTTTTTATTGCTGCAAGTCTGCTTGCAGCTGGCTTTCCAGCAGTCATTAATGATTTACAGCTACATCTTGAGGTACAGCAGTCGTTGAGGTACACGCAGTGGTCATATGTGTGCTTTGTGCAACTGTGTTTCTGTTATGTATTTGCTTTTGGTTTACCCTTGGCAGCTATGTAGCTAAGCTGCATGTAGGAGTCACCTTCATCGTTGAGGGACAGCTGGGTTTTGAAGCTGTGGGAATGGTAGTTGGATAGTGTGTGACAAAACCTTACAGATGCAAAGCAACAGAGCTAATTTTGTACAGTTTGTTATCATTGCAGCAGGAGGATTTGGAATAGTTGTTAGAACAGTTAGCTATAAAACCTATAGGTACCAGTTAATATTTGGGATTTGTGGCATTTGAGCCAGCAATAGGAGCTATGTGCTAACCCAGCTGTTGTTAGGTTTATTGTAAAGATGTGCTGGTATTGATGGGAATGTAGGAGGTGCACAATGCTTACAATCATCATAACAAAGCATACGAGATGTAGCTAAATGTAGATGTTTGCTTTTGGGTTCAGGGGTACATGTGAAGGAAATATGTGTCACATAGGTGGTTGTGTAGATATATTGATGCACCCATAGAAATGGAATGCAGGAACTGGACCCAGTGTAGTCATGGTTTTAATAATAGCCAGGCAGAGATGGAGGCAGCAACAGCATCTGTTGCTGATGTAAAGAAAAGGAAGTCATTAAGAGTTTGTTTGAGATCATGAAACAGCTAACACAGCAGAGGGCAAAGCCAGAGAGGTA',
                        'primers': {'forward': 'ATGCAACAGTTATGCACAGA', 'reverse': 'ACTGCGTGTACCTCAACGAC'},
                        'probe': {'sequence': 'CAGCTGCAAACAACTGCAATGT', 'dye': 'HEX', 'quencher': 'BHQ1'}
                    }
                }
            }

    def load_fasta(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Open FASTA File", "", "FASTA Files (*.fasta *.fa *.txt)")
        if file_name:
            with open(file_name) as file:
                content = file.readlines()
                sequence = "".join([line.strip() for line in content if not line.startswith(">")])
                self.dna_input.setPlainText(sequence)

    def load_example(self):
        self.dna_input.setPlainText(self.db['targets']['SARS-CoV-2']['sequence'])
        self.forward_primer.setText(self.db['targets']['SARS-CoV-2']['primers']['forward'])
        self.reverse_primer.setText(self.db['targets']['SARS-CoV-2']['primers']['reverse'])
        self.probe_sequence.setText(self.db['targets']['SARS-CoV-2']['probe']['sequence'])
        self.dye_combo.setCurrentText(self.db['targets']['SARS-CoV-2']['probe']['dye'])

    def clear_fields(self):
        self.dna_input.clear()
        self.forward_primer.clear()
        self.reverse_primer.clear()
        self.probe_sequence.clear()
        self.dye_combo.setCurrentIndex(0)
        self.multiplex_check.setChecked(False)
        self.cycles_input.setText("40")
        self.efficiency_input.setText("0.95")
        self.threshold_input.setText("0.1")
        self.denature_temp.setText("95")
        self.anneal_temp.setText("60")
        self.extend_temp.setText("72")
        self.ax.clear()
        self.canvas.draw()
        self.results_text.clear()

    def validate_inputs(self):
        try:
            cycles = int(self.cycles_input.text())
            efficiency = float(self.efficiency_input.text())
            threshold = float(self.threshold_input.text())
            
            if not (0 < efficiency <= 1):
                raise ValueError("Efficiency must be between 0 and 1")
                
            if cycles <= 0:
                raise ValueError("Cycle count must be positive")
                
            if threshold <= 0:
                raise ValueError("Threshold must be positive")
                
            if not self.dna_input.toPlainText().strip():
                raise ValueError("DNA sequence is required")
                
            if not self.forward_primer.text().strip() or not self.reverse_primer.text().strip():
                raise ValueError("Both forward and reverse primers are required")
                
            return True
            
        except ValueError as e:
            QMessageBox.warning(self, "Input Error", str(e))
            return False

    def run_simulation(self):
        if not self.validate_inputs():
            return
            
        try:
            cycles = int(self.cycles_input.text())
            efficiency = float(self.efficiency_input.text())
            threshold = float(self.threshold_input.text())
            
            # Simulate PCR amplification
            fluorescence = []
            product = 1
            ct_value = None
            
            for i in range(cycles):
                product *= (1 + efficiency)
                noise = random.uniform(-0.1, 0.1)
                current_fluo = product + noise
                fluorescence.append(current_fluo)
                
                # Determine Ct value (first cycle where fluorescence crosses threshold)
                if ct_value is None and current_fluo > threshold:
                    ct_value = i
            
            # Plot results
            self.ax.clear()
            self.ax.plot(range(cycles), fluorescence, color='#1f77b4', linewidth=2)
            self.ax.axhline(y=threshold, color='r', linestyle='--', label='Threshold')
            
            if ct_value is not None:
                self.ax.axvline(x=ct_value, color='g', linestyle=':', label=f'Ct = {ct_value}')
            
            self.ax.set_title("PCR Amplification Curve", fontsize=12)
            self.ax.set_xlabel("Cycle Number", fontsize=10)
            self.ax.set_ylabel("Fluorescence (RFU)", fontsize=10)
            self.ax.grid(True, linestyle='--', alpha=0.7)
            self.ax.legend()
            self.canvas.draw()
            
            # Display results text
            result_text = f"""=== PCR Simulation Results ===
Cycles: {cycles}
Efficiency: {efficiency:.2f}
Threshold: {threshold:.2f}
Ct value: {ct_value if ct_value is not None else "Not reached"}

=== Thermal Profile ===
Denaturation: {self.denature_temp.text()}°C
Annealing: {self.anneal_temp.text()}°C
Extension: {self.extend_temp.text()}°C

=== Primer/Probe Info ===
Forward: {self.forward_primer.text()}
Reverse: {self.reverse_primer.text()}
Probe: {self.probe_sequence.text()} ({self.dye_combo.currentText()})
Multiplex: {'Yes' if self.multiplex_check.isChecked() else 'No'}"""
            
            self.results_text.setPlainText(result_text)
            
        except Exception as e:
            QMessageBox.critical(self, "Simulation Error", f"An error occurred: {str(e)}")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = PCRSimulator()
    window.show()
    sys.exit(app.exec_())