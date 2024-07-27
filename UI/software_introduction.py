from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QResizeEvent
from PySide6.QtWidgets import QWidget, QVBoxLayout, QSizePolicy, QSpacerItem, QHBoxLayout

from qfluentwidgets import TitleLabel, BodyLabel, HyperlinkLabel
from qfluentwidgets import TextWrap, getFont


class IntroducationInterface(QWidget):
    def __init__(self, object_name: str, parent: QWidget = None) -> None:
        super().__init__(parent)
        self._init_ui()
        self.setStyleSheet(
            """
            QWidget{
                background-color: #fafafa;
                border-radius: 10px;
                margin: 10px;
            }
            """
        )
        self.setObjectName(object_name.replace(' ', '-'))

    def _init_ui(self):
        self.content = 'Repeat-Enhancing Featured Ion-Guided Stoichiometry (RE-FIGS) is a complete and compact solution on DI-SPA data for more confident identifications and corresponding label free quantifications.'
        self.introduction_title_label = TitleLabel(
            'Simple Introduction', self)
        self.cite_label = TitleLabel(
            'Cite Paper', self
        )
        self.content_label = BodyLabel(self)
        self.content_label.setSizePolicy(
            QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        self.content_label.setFont(getFont(fontSize=14))

        self.request_label = BodyLabel(self)
        self.request_label.setText(
            'If this software useful to you, please cite our paper!')
        self.paper_url = HyperlinkLabel(
            QUrl('https://pubs.acs.org/doi/10.1021/acs.jproteome.3c00050'), 'paper')

        self.vboxlayout = QVBoxLayout()
        self.hboxlayout = QHBoxLayout()
        self.setLayout(self.vboxlayout)
        self.vboxlayout.setAlignment(Qt.AlignmentFlag.AlignLeft)
        self.vboxlayout.addWidget(self.introduction_title_label)
        self.vboxlayout.addWidget(self.content_label)
        self.vboxlayout.addWidget(self.cite_label)
        self.vboxlayout.addLayout(self.hboxlayout)
        self.vboxlayout.setContentsMargins(0, 10, 0, 10)
        self.vboxlayout.addItem(QSpacerItem(
            2, 10, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding))

        self.hboxlayout.addWidget(self.request_label)
        self.hboxlayout.addWidget(self.paper_url)
        self.hboxlayout.setAlignment(Qt.AlignmentFlag.AlignLeft)
        self.hboxlayout.setContentsMargins(0, 0, 5, 0)

    def resizeEvent(self, event: QResizeEvent) -> None:
        self._adjustText()

    def _adjustText(self):
        w = self.width()
        chars = max(min(w / 8, 200), 30)

        self.content_label.setText(
            TextWrap.wrap(self.content, chars, False)[0])
