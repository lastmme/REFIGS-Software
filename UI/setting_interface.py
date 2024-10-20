from typing import Union

from qfluentwidgets import (
    SettingCard,
    SpinBox,
    DoubleSpinBox,
    OptionsValidator,
    ComboBoxSettingCard,
    OptionsConfigItem,
    TitleLabel,
    ConfigItem,
    FluentIconBase,
    SettingCardGroup,
    ExpandSettingCard,
    PrimaryPushButton,
    PushButton,
    SwitchSettingCard,
    Path,
    Dialog,
    RangeConfigItem,
    RangeValidator,
    RangeSettingCard,
    PushSettingCard,
    SwitchSettingCard,
    ScrollArea,
    VBoxLayout,
    MessageBox
)
from qfluentwidgets import FluentIcon as FIF
from qfluentwidgets.components.settings.folder_list_setting_card import FolderItem
from PySide6.QtWidgets import QWidget, QFileDialog, QSizePolicy, QHBoxLayout,QMessageBox
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QIcon

import pixmaps
from multiprocessing import cpu_count


class DoubleSpinCard(SettingCard):
    def __init__(
        self,
        icon: Union[str, QIcon, FluentIconBase],
        title: str,
        content: str,
        min_value: float,
        max_value: float,
        step_value: float,
        parent=None
    ):
        super().__init__(icon, title, content, parent)
        # spin_box
        self.spin_box = DoubleSpinBox(self)
        self.spin_box.setFixedWidth(150)
        self.spin_box.setMinimum(min_value)
        self.spin_box.setMaximum(max_value)
        self.spin_box.setSingleStep(step_value)
        # button
        self.hBoxLayout.addWidget(
            self.spin_box, 0, Qt.AlignmentFlag.AlignRight)
        self.hBoxLayout.addSpacing(16)

    def set_content(self, content: str):
        self.content_label.setText(content)

    @property
    def value(self):
        return self.spin_box.value()


class SpinCard(SettingCard):
    def __init__(
        self,
        icon: Union[str, QIcon, FluentIconBase],
        title: str,
        content: str,
        min_value: float,
        max_value: float,
        step_value: float,
        parent=None
    ):
        super().__init__(icon, title, content, parent)
        # spin_box
        self.spin_box = SpinBox(self)
        self.spin_box.setFixedWidth(150)
        self.spin_box.setMinimum(min_value)
        self.spin_box.setMaximum(max_value)
        self.spin_box.setSingleStep(step_value)
        # button
        self.hBoxLayout.addWidget(
            self.spin_box, 0, Qt.AlignmentFlag.AlignRight)
        self.hBoxLayout.addSpacing(16)

    def set_content(self, content: str):
        self.content_label.setText(content)
    def setMaxvalue(self, files):
        value=min(cpu_count(),len(files))
        self.spin_box.setMaximum(value)

    @property
    def value(self):
        return self.spin_box.value()


class FileListSettingCard(ExpandSettingCard):
    """ File list setting card """

    fileChanged = Signal(list)

    def __init__(self, icon, title: str, filter: str, content: str = None, directory="./", parent=None):
        super().__init__(icon, title, content, parent)
        self._dialogDirectory = directory
        self._filter = filter
        self.addFolderButton = PushButton(
            self.tr('添加文件'), self, FIF.FOLDER_ADD
        )
        self.deleteAllFilesButton = PushButton(
            self.tr('删除全部文件'), self, FIF.DELETE
        )

        self.files = []
        self.items = []
        self.__initWidget()

    def __initWidget(self):
        self.addWidget(self.addFolderButton)
        self.addWidget(self.deleteAllFilesButton)

        # initialize layout
        self.viewLayout.setSpacing(0)
        self.viewLayout.setAlignment(Qt.AlignTop)
        self.viewLayout.setContentsMargins(0, 0, 0, 0)
        for file in self.files:
            self.__addFolderItem(file)

        self.deleteAllFilesButton.clicked.connect(self.__showDeleteDialog)
        self.addFolderButton.clicked.connect(self.__showFolderDialog)

    def __showDeleteDialog(self):
        title = self.tr('Are you sure you want to delete all files?')
        w = Dialog(title, '', self.window())
        w.yesSignal.connect(self.deleteAllFiles)
        w.exec_()

    def deleteAllFiles(self):
        """ Delete all files """
        for i in reversed(range(self.viewLayout.count())):
            item = self.viewLayout.itemAt(i).widget()
            if isinstance(item, FolderItem):
                self.__removeFolder(item)

    def __showFolderDialog(self):
        """ show folder dialog """
        files, _ = QFileDialog.getOpenFileUrls(
            self, self.tr("选择需要分析的数据"), self._dialogDirectory, self._filter
        )

        for file in files:
            path = file.toLocalFile()
            if path not in self.files:
                self.__addFolderItem(path)
                self.files.append(path)
                self.fileChanged.emit(self.files)

    def __addFolderItem(self, folder: str):
        """ add file item """
        item = FolderItem(folder, self.view)
        item.removed.connect(self.__showConfirmDialog)
        self.viewLayout.addWidget(item)
        item.show()
        self._adjustViewSize()

    def __showConfirmDialog(self, item: FolderItem):
        """ show confirm dialog """
        name = Path(item.folder).name
        title = self.tr('Are you sure you want to delete the file?')
        content = self.tr("If you delete the ") + f'"{name}"' + \
            self.tr(" file and remove it from the list, the file will no "
                    "longer appear in the list, but will not be deleted.")
        w = Dialog(title, content, self.window())
        w.yesSignal.connect(lambda: self.__removeFolder(item))
        w.exec_()

    def __removeFolder(self, item: FolderItem):
        """ remove folder """
        if item.folder not in self.files:
            return

        self.files.remove(item.folder)
        self.viewLayout.removeWidget(item)
        item.deleteLater()
        self._adjustViewSize()

        self.fileChanged.emit(self.files)


class PointSliderCard(RangeSettingCard):
    def __init__(self, configItem: RangeConfigItem, icon: Union[str, QIcon, FluentIconBase], title, content=None, parent=None):
        super().__init__(configItem, icon, title, content, parent)
        self.max_value = configItem.validator.max
        self.true_value = configItem.value / self.max_value
        self.valueLabel.setNum(self.true_value)

    def setValue(self, value):
        super().setValue(value)
        self.true_value = value / self.max_value
        self.valueLabel.setNum(self.true_value)
        self.valueLabel.adjustSize()


class SettingInterface(ScrollArea):

    excute_signal = Signal()

    def __init__(self, object_name: str, parent=None):
        super().__init__(parent)
        self.setObjectName(object_name)
        self.scroll_widget = QWidget()
        self.expand_layout = VBoxLayout(self.scroll_widget)

        # Setting label
        self.setting_label = TitleLabel(self.tr("Settings"), self)

        # spectrum data and library position
        self.soruce_groups = SettingCardGroup(
            self.tr('原始数据导入'),
            self.scroll_widget
        )
        self.result_file_card = PushSettingCard(
            self.tr('选择文件夹'),
            FIF.FOLDER,
            '最终结果的存储位置',
            '',
            self.soruce_groups
        )
        self.result_file_card.contentLabel.text()
        self.library_file_card = PushSettingCard(
            self.tr('选择文件'),
            ':/normal/svgs/normal/library.svg',
            '图谱库文件路径',
            '',
            self.soruce_groups
        )
        self.mzxml_files_card = FileListSettingCard(
            icon=':/normal/svgs/normal/mzxml.svg',
            title=self.tr('选择需要分析的 mzXML 文件'),
            directory='./',
            filter='*.mzXML',
            parent=self.soruce_groups
        )
        self.cpus=cpu_count()
        self.worker_num_card = SpinCard(
            icon=':/normal/svgs/normal/amount.svg',
            title='选择 CPU 核心数',
            content='',
            min_value=1,
            max_value=self.cpus,
            step_value=1,
            parent=self.soruce_groups
        )
        self.mzxml_files_card.fileChanged.connect(self.worker_num_card.setMaxvalue)
        # self.mzml_files_card = FileListSettingCard(
        #     icon=':/normal/svgs/normal/mzml.svg',
        #     title=self.tr('选择需要分析的 mzML 文件'),
        #     directory='./',
        #     filter='*.mzML',
        #     parent=self.soruce_groups
        # )
        self.mzxml_files_card.addFolderButton.setText(self.tr('添加文件'))
        # self.mzml_files_card.addFolderButton.setText(self.tr('添加文件'))

        # decoy library parameters
        self.decoy_parameters_groups = SettingCardGroup(
            self.tr('诱饵库参数设置'),
            self.scroll_widget
        )
        self.strategy_card = ComboBoxSettingCard(
            OptionsConfigItem(
                'decoy_parameters',
                'strategy',
                '陷阱-诱饵策略',
                OptionsValidator(
                    [
                        True,
                        False
                    ]
                ),
            ),
            ':/normal/svgs/normal/strategy.svg',
            self.tr('诱饵库策略选择'),
            '',
            [
                self.tr('trap-decoy'),
                self.tr('decoy')
            ],
            self.decoy_parameters_groups
        )
        self.entrap_distance_card = RangeSettingCard(
            RangeConfigItem(
                'decoy_parameters',
                'entrapment_distance',
                23,
                RangeValidator(10, 100),
            ),
            ':/normal/svgs/normal/distance.svg',
            self.tr('entrapment_distance'),
            parent=self.decoy_parameters_groups
        )
        self.decoy_distance_card = RangeSettingCard(
            RangeConfigItem(
                'decoy_parameters',
                'decoy_distance',
                33,
                RangeValidator(10, 100),
            ),
            ':/normal/svgs/normal/distance.svg',
            self.tr('decoy_distance'),
            parent=self.decoy_parameters_groups
        )
        self.fdr_card = DoubleSpinCard(
            icon=FIF.CARE_RIGHT_SOLID,
            title='fdr 阈值设置',
            content='',
            min_value=0,
            max_value=1,
            step_value=0.01,
            parent=self.decoy_parameters_groups
        )
        # Csodiaq Parameters
        self.csodiaq_paramter_groups = SettingCardGroup(
            self.tr('Csodiaq 参数设置'),
            self.scroll_widget
        )
        self.tolerance_card = RangeSettingCard(
            RangeConfigItem(
                'csodiaq_parameters',
                'tolerance',
                20,
                RangeValidator(0, 50),
            ),
            FIF.SEARCH,
            self.tr('tolerance'),
            parent=self.csodiaq_paramter_groups
        )
        self.correted_switch_card = SwitchSettingCard(
            ':/normal/svgs/normal/switch.svg',
            self.tr('is_corrected'),
            '是否采用矫正',
            ConfigItem(
                'csodiaq_parameters',
                'is_correct',
                True
            ),
            self.csodiaq_paramter_groups
        )
        # self.spectra_pool_card = RangeSettingCard(
        #     RangeConfigItem(
        #         'csodiaq_parameters',
        #         'max_num_spectra_pool',
        #         30,
        #         RangeValidator(0, 100),
        #     ),
        #     ':/normal/svgs/normal/amount.svg',
        #     self.tr('maxQuerySpectraToPool'),
        #     parent=self.csodiaq_paramter_groups
        # )
        # RE-FIGS 定性参数设置
        self.refigs_identify_groups = SettingCardGroup(
            self.tr('RE-FIGS 参数设置'),
            self.scroll_widget
        )
        self.top10_card = SwitchSettingCard(
            ':/normal/svgs/normal/switch.svg',
            self.tr('是否采用 Top 10 的峰'),
            '',
            ConfigItem(
                'RE-FIGS_identification_parameters',
                'is_top_10',
                True
            ),
            self.csodiaq_paramter_groups
        )
        self.start_cycle_card = RangeSettingCard(
            RangeConfigItem(
                'RE-FIGS_identification_parameters',
                'start_cycle',
                1,
                RangeValidator(1, 10),
            ),
            ':/normal/svgs/normal/start_cycle.svg',
            self.tr('start_cycle'),
            self.tr('the cycle No. starts from.'),
            self.csodiaq_paramter_groups
        )
        self.end_cycle_card = RangeSettingCard(
            RangeConfigItem(
                'RE-FIGS_identification_parameters',
                'end_cycle',
                14,
                RangeValidator(10, 20),
            ),
            ':/normal/svgs/normal/end_cycle.svg',
            self.tr('end_cycle'),
            self.tr('the cycle No. ends with'),
            self.csodiaq_paramter_groups
        )
        
        self.auto_search_refigs_params=SwitchSettingCard(
            ':/normal/svgs/normal/switch.svg',
            self.tr('自动搜索 RE-FIGS 参数'),
            '是否自动搜索 RE-FIGS 参数：peaks_shared_limit, cosine score limit, sqrt cosine score, cycle within limit',
            ConfigItem(
                'RE-FIGS_identification_parameters',
                'auto_search_refigs_params',
                True
            ),
            self.csodiaq_paramter_groups
        )
        
        self.good_shared_limit_card = RangeSettingCard(
            RangeConfigItem(
                'RE-FIGS_identification_parameters',
                'good_shared_limit',
                5,
                RangeValidator(3, 10),
            ),
            ':/normal/svgs/normal/shared_count.svg',
            self.tr('peaks_shared_limit'),
            self.tr('the threshold to select good target'),
            self.csodiaq_paramter_groups
        )
        self.good_cos_sim_limit_card = PointSliderCard(
            RangeConfigItem(
                'RE-FIGS_identification_parameters',
                'good_cos_sim_limit',
                80,
                RangeValidator(50, 100),
            ),
            ':/normal/svgs/normal/cosine.svg',
            self.tr('cosine score limit'),
            self.tr('the threshold to select good target. range(0,1)'),
            self.csodiaq_paramter_groups
        )
        self.good_cos_sim_limit_card.slider.setSingleStep(5)

        self.good_sqrt_cos_sim_limit_card = PointSliderCard(
            RangeConfigItem(
                'RE-FIGS_identification_parameters',
                'good_sqrt_cos_sim_limit',
                90,
                RangeValidator(50, 100),
            ),
            ':/normal/svgs/normal/sqrt_cosine.svg',
            self.tr('sqrt cosine score'),
            self.tr('the threshold to select good target. range(0,1)'),
            self.csodiaq_paramter_groups
        )
        self.good_sqrt_cos_sim_limit_card.slider.setSingleStep(5)

        self.good_count_within_cycle_limit = RangeSettingCard(
            RangeConfigItem(
                'RE-FIGS_identification_parameters',
                'good_count_within_cycle_limit',
                5,
                RangeValidator(1, 10),
            ),
            ':/normal/svgs/normal/amount.svg',
            self.tr('cycle within limit'),
            self.tr('the threshold to select good target'),
            parent=self.csodiaq_paramter_groups
        )
        # self.scans_per_cycle = RangeSettingCard(
        #     RangeConfigItem(
        #         'RE-FIGS_identification_parameters',
        #         'scans_per_cycle',
        #         1500,
        #         RangeValidator(500, 3000),
        #     ),
        #     ':/normal/svgs/normal/amount.svg',
        #     self.tr('scan per cycle'),
        #     self.tr('the scan number in each cycle.'),
        # )
        # self.scans_per_cycle.slider.setSingleStep(100)

        self.seed_card = RangeSettingCard(
            RangeConfigItem(
                'RE-FIGS_identification_parameters',
                'seed',
                13,
                RangeValidator(0, 100),
            ),
            ':/normal/svgs/normal/random_seed.svg',
            self.tr('random seed'),
            self.tr('the seed to randomly choose decoy'),
            self.csodiaq_paramter_groups
        )
        # excute button
        self.qhboxlayout = QHBoxLayout()

        self.excute_button = PrimaryPushButton(
            self.tr('开始分析'), self.scroll_widget)
        self.excute_button.setMinimumHeight(30)
        self.excute_button.setMaximumWidth(500)
        self.excute_button.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Minimum
        )

        self._init_ui()
        self._set_qss()
        self._connect_signal2solt()

    def _init_ui(self):
        self.setting_label.move(60, 63)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setViewportMargins(0, 120, 0, 20)
        self.setWidget(self.scroll_widget)
        self.setWidgetResizable(True)
        # source_groups
        self.soruce_groups.addSettingCards(
            [
                self.result_file_card,
                self.library_file_card,
                self.mzxml_files_card,
                self.worker_num_card
                # self.mzml_files_card
            ]
        )
        self.decoy_parameters_groups.addSettingCards(
            [
                self.strategy_card,
                self.entrap_distance_card,
                self.decoy_distance_card,
                self.fdr_card
            ]
        )
        # csodiaq_paramter_groups
        self.csodiaq_paramter_groups.addSettingCards(
            [
                self.tolerance_card,
                self.correted_switch_card,
                # self.spectra_pool_card
            ]
        )
        # RE-FIGS Identification
        self.refigs_identify_groups.addSettingCards(
            [
                self.top10_card,
                self.start_cycle_card,
                self.end_cycle_card,
                self.auto_search_refigs_params,
                self.good_shared_limit_card,
                self.good_cos_sim_limit_card,
                self.good_sqrt_cos_sim_limit_card,
                self.good_count_within_cycle_limit,
                # self.scans_per_cycle,
                self.seed_card
            ]
        )

        # expand layout
        self.expand_layout.setSpacing(28)
        self.expand_layout.setContentsMargins(60, 10, 60, 0)
        self.expand_layout.addWidget(self.soruce_groups)
        self.expand_layout.addWidget(self.decoy_parameters_groups)
        self.expand_layout.addWidget(self.csodiaq_paramter_groups)
        self.expand_layout.addWidget(self.refigs_identify_groups)

        # qhboxlayout
        self.qhboxlayout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.qhboxlayout.addWidget(self.excute_button)
        self.qhboxlayout.setContentsMargins(0, 2, 10, 2)

        self.expand_layout.addLayout(self.qhboxlayout)

    def _set_qss(self):
        self.scroll_widget.setObjectName('scrollWidget')
        self.setting_label.setObjectName('settingLabel')
        self.setStyleSheet(
            """
            QScrollArea {
                background-color: rgb(249, 249, 249);
                border-radius: 10px;
                border: none;
            }
            """)
        self.scroll_widget.setStyleSheet(
            """
            QWidget#scrollWidget {
                background-color: rgb(249, 249, 249);
                border-radius: 10px;

            }
            """
        )
        self.setting_label.setStyleSheet(
            """
            QLabel#settingLabel {
                background-color: transparent;
            }
            """
        )

    def _connect_signal2solt(self):
        self.library_file_card.clicked.connect(self._choose_library)
        self.result_file_card.clicked.connect(self._choose_reuslt_folder)
        self.excute_button.clicked.connect(self.excute_signal)

    def _choose_reuslt_folder(self):
        folder = QFileDialog.getExistingDirectoryUrl(
            self,
            self.tr('选择结果文件夹'),
            './',
        ).toLocalFile()
        self.result_file_card.setContent(folder)

    def _choose_library(self):
        file, _ = QFileDialog.getOpenFileUrl(
            self,
            self.tr('选择图谱库文件'),
            './',
            '*.mgf'
        )
        self.library_file_card.setContent(file.toLocalFile())

    def reset(self):
        """
            TODO:  将所有设置恢复为默认值
        """
        self.mzxml_files_card.deleteAllFiles()
        self.result_file_card.contentLabel.setText('')
        self.result_file_card.contentLabel.setVisible(bool(''))
        self.library_file_card.contentLabel.setText('')
        self.library_file_card.contentLabel.setVisible(bool('')) 
        self.worker_num_card.spin_box.setMaximum(5)
        self.worker_num_card.spin_box.setValue(1)
        self.strategy_card.configItem.value = True
        self.entrap_distance_card.configItem.value = 23
        self.decoy_distance_card.configItem.value = 33
        self.fdr_card.spin_box.setValue(0)
        self.tolerance_card.configItem.value = 20
        self.correted_switch_card.configItem.value = True
        # self.spectra_pool_card.configItem.value = 30
        self.top10_card.configItem.value = True
        self.start_cycle_card.configItem.value = 1
        self.end_cycle_card.configItem.value = 14
        self.auto_search_refigs_params.configItem.value = True
        self.good_shared_limit_card.configItem.value = 5
        self.good_cos_sim_limit_card.configItem.value = 80
        self.good_sqrt_cos_sim_limit_card.configItem.value = 90
        self.good_count_within_cycle_limit.configItem.value = 5
        # self.scans_per_cycle.configItem.value = 1500
        self.seed_card.configItem.value = 13
        self.excute_button.setEnabled(True)

    def get_params(self):
        self.excute_button.setEnabled(False)
        if len(self.mzxml_files_card.files) == 0:
            self.excute_button.setEnabled(True)
            self.show_error_message('Error', 'Please select mzXML files')
            return None
        if self.library_file_card.contentLabel.text() == '':
            self.excute_button.setEnabled(True)
            self.show_error_message('Error', 'Please select library file')
            return None
        if self.result_file_card.contentLabel.text() == '':
            self.excute_button.setEnabled(True)
            self.show_error_message('Error', 'Please select result folder')
            return None
        
        return self.worker_num_card.value, {
            'mzxml_files': self.mzxml_files_card.files,
            # 'mzml_files': self.mzml_files_card.files,
            'result_dir': self.result_file_card.contentLabel.text(),
            'library_file': self.library_file_card.contentLabel.text(),
            'is_entrap': self.strategy_card.configItem.value,
            'fdr_value': self.fdr_card.value,
            'entrap_distance': self.entrap_distance_card.configItem.value,
            'decoy_distance': self.decoy_distance_card.configItem.value,
            'tolerance': self.tolerance_card.configItem.value,
            'is_corrected': self.correted_switch_card.configItem.value,
            'is_top10': self.top10_card.configItem.value,
            'start_cycle': self.start_cycle_card.configItem.value,
            'end_cycle': self.end_cycle_card.configItem.value,
            'is_auto_search_refigs_params': self.auto_search_refigs_params.configItem.value,
            'good_shared_limit': self.good_shared_limit_card.configItem.value,
            'good_cos_sim_limit': float(self.good_cos_sim_limit_card.configItem.value) / 100,
            'good_sqrt_cos_sim_limit': float(self.good_sqrt_cos_sim_limit_card.configItem.value)/100,
            'good_count_within_cycle_limit': self.good_count_within_cycle_limit.configItem.value,
            # 'scans_per_cycle': self.scans_per_cycle.configItem.value,
            'seed': self.seed_card.configItem.value
        }
    
    def show_error_message(self, title, message):
        msg_box = QMessageBox()
        msg_box.setIcon(QMessageBox.Critical)
        msg_box.setWindowTitle(title)
        msg_box.setText(message)
        msg_box.setStandardButtons(QMessageBox.Ok)
        msg_box.exec_()
