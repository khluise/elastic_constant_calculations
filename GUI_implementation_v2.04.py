import sys
import os
import subprocess
from PyQt5.QtWidgets import QDialog, QApplication, QWidget, QPushButton, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QFileDialog, QMessageBox, QTextBrowser
from PyQt5.QtGui import QFont, QDesktopServices, QIcon
from PyQt5.QtCore import Qt, QUrl

class MyWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle("Elastic Calculation")
        self.setGeometry(100, 100, 600, 250)

        layout = QVBoxLayout()

        self.castep_input = self.create_file_input("Input the .castep file", "*.castep")
        layout.addLayout(self.castep_input['layout'])

        self.elastic_input = self.create_file_input("Input the elastic constant.txt file", "*.txt")
        layout.addLayout(self.elastic_input['layout'])

        self.compound_input = self.create_text_input("Input the compound")
        layout.addLayout(self.compound_input['layout'])

        self.output_input = self.create_file_input("Output file name", "*.txt", save=True)
        layout.addLayout(self.output_input['layout'])

        button_layout = QHBoxLayout()

        about_button = QPushButton("About me")
        about_button.clicked.connect(self.show_about)
        about_button.setFont(QFont("Arial", 12))
        button_layout.addWidget(about_button)

        calculate_button = QPushButton("Calculate")
        calculate_button.setFont(QFont("Arial", 12))
        calculate_button.clicked.connect(self.calculate)
        button_layout.addWidget(calculate_button)

        credits_button = QPushButton("Credits")
        credits_button.clicked.connect(self.show_credits)
        credits_button.setFont(QFont("Arial", 12))
        button_layout.addWidget(credits_button)

        layout.addLayout(button_layout)
        self.setLayout(layout)

    def calculate(self):
        castep_file_path = self.castep_input['line_edit'].text()
        elastic_file_path = self.elastic_input['line_edit'].text()
        compound_name_formula = self.compound_input['line_edit'].text()
        output_file_path = self.output_input['line_edit'].text()
        if not castep_file_path or not elastic_file_path or not compound_name_formula or not output_file_path:
            QMessageBox.warning(self, "Input Error", "Please select all input files and specify the output file.")
            return
        self.write_output_to_file(castep_file_path, elastic_file_path, compound_name_formula, output_file_path)
    
    def write_output_to_file(self, castep_file_path, elastic_file_path, compound_name_formula, output_file_path):
        with open(output_file_path, 'w', encoding="utf-8") as file:
            sys.stdout = file
            import elastic_calculation
            elastic_calculation.calculate(castep_file_path, elastic_file_path, compound_name_formula)
            file.close()
        self.show_completion_message(output_file_path)
        sys.stdout = sys.__stdout__
        
    def show_completion_message(self, file_path):
        msg = QMessageBox(self)
        msg.setWindowTitle("Success")
        msg.setText("Calculation saved successfully!")
        open_button = msg.addButton("Open File", QMessageBox.AcceptRole)
        back_button = msg.addButton("Go Back", QMessageBox.RejectRole)

        msg.exec_()

        if msg.clickedButton() == open_button:
            self.open_file(file_path)
        elif msg.clickedButton() == back_button:
            self.close()

    def open_file(self, file_path):
        """ Opens the saved file in the default system application """
        try:
            if sys.platform.startswith("win"):
                os.startfile(file_path)
            elif sys.platform.startswith("darwin"):
                subprocess.call(["open", file_path])
            else:
                subprocess.call(["xdg-open", file_path])
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Could not open file:\n{e}")


    def create_file_input(self, label_text, file_filter, save=False):
        layout = QHBoxLayout()
        label = QLabel(label_text)
        label.setFont(QFont("Arial", 12))
        layout.addWidget(label)
        line_edit = QLineEdit()
        line_edit.setFont(QFont("Arial", 12))
        layout.addWidget(line_edit)
        browse_button = QPushButton("Browse file")
        browse_button.setFont(QFont("Arial", 12))
        browse_button.clicked.connect(lambda: self.browse_file(line_edit, file_filter, save))
        layout.addWidget(browse_button)
        return {'layout': layout, 'line_edit': line_edit}
    
    def create_text_input(self, label_text):
        layout = QHBoxLayout()
        label = QLabel(label_text)
        label.setFont(QFont("Arial", 12))
        layout.addWidget(label)
        line_edit = QLineEdit()
        line_edit.setFont(QFont("Arial", 12))
        layout.addWidget(line_edit)
        return {'layout': layout, 'line_edit': line_edit}

    def browse_file(self, line_edit, file_filter, save):
        if save:
            file_name, _ = QFileDialog.getSaveFileName(self, "Select output file", filter=file_filter)
        else:
            file_name, _ = QFileDialog.getOpenFileName(self, "Select file", filter=file_filter)
        if file_name:
            line_edit.setText(file_name)

    def show_about(self):
        about_text = """
        <p>This application was developed to process and analyze different
        elastic and thermophysical properties using the elastic constants
        .txt file generated from the CASTEP software.</p>

        <p>For the formulas used in the calculations, please refer to the
        credit section. If you want to add more calculation or find any
        formulas to be problematic please contact with me.</p>
        <p>
        <i>Developer:</i> <b>Md. Kamrul Hassan Luise</b><br>
        <i>Contact:</i> <a href='mailto:luisekh076@gmail.com'>luisekh076@gmail.com</a><br>
        <i>Phone no:</i> +8801315188626
        </p>
        """
        
        dialog = QDialog(self)
        dialog.setWindowTitle("About Me")
        dialog.setWindowIcon(QIcon("about_me.ico"))
        dialog.setMinimumSize(500, 300)

        layout = QVBoxLayout(dialog)
        
        aboutme_box = QTextBrowser(dialog)
        aboutme_box.setHtml(about_text)
        aboutme_box.setReadOnly(True)
        aboutme_box.setTextInteractionFlags(Qt.TextSelectableByMouse)
        aboutme_box.setOpenExternalLinks(True) 
        aboutme_box.anchorClicked.connect(self.open_email_client) 
        layout.addWidget(aboutme_box)
        
        font = QFont("Arial", 10)
        aboutme_box.setFont(font)

        close_button = QPushButton("Close", dialog)
        close_button.setFont(QFont("Arial", 12))
        close_button.clicked.connect(dialog.accept)
        layout.addWidget(close_button)
        
        dialog.setLayout(layout)
        dialog.exec_()

    def open_email_client(self, url: QUrl):
        # Open the default email client when the mailto link is clicked
        if url.scheme() == "mailto":
            QDesktopServices.openUrl(url)

    def show_credits(self):
        credit_text = """
        <P>This project was done with the help of:</P>
        <p>
        <b>Jubair Hossain Abir</b><br>
        <i>Contact:</i> <a href='mailto:abir100r8@gmail.com'>abir100r8@gmail.com</a><br>
        <i>Phone no:</i> +8801866457810
        </p>

        <p>If you want to check the folmulas used in this program <br>
        click the button below.</p>
        """
        dialog = QDialog(self)
        dialog.setWindowTitle("Credits")
        dialog.setWindowIcon(QIcon("credits.ico"))  # Optional: Add an icon
        dialog.setMinimumSize(500, 50)

        layout = QVBoxLayout(dialog)
        
        credits_box = QTextBrowser(dialog)
        credits_box.setHtml(credit_text)
        credits_box.setReadOnly(True)
        credits_box.setTextInteractionFlags(Qt.TextSelectableByMouse)
        credits_box.setOpenExternalLinks(True) 
        credits_box.anchorClicked.connect(self.open_email_client) 
        layout.addWidget(credits_box)
        
        font = QFont("Arial", 10)
        credits_box.setFont(font)

        check_formulas_button = QPushButton("Check Formulas")
        check_formulas_button.setFont(QFont("Arial", 12))
        check_formulas_button.clicked.connect(self.open_docx)
        layout.addWidget(check_formulas_button)

        dialog.setLayout(layout)
        dialog.exec_()

    def open_docx(self):
        docx_path = "Elastic constants calculation formula.docx"  # Replace with the actual path to the DOCX file
        if os.path.exists(docx_path):
            os.startfile(docx_path) if sys.platform == "win32" else os.system(f"open {docx_path}")
        else:
            QMessageBox.warning(self, "Error", "File not found.")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MyWindow()
    window.show()
    sys.exit(app.exec_())
