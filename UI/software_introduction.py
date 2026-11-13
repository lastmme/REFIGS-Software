"""In-app user guide. Content describes user operations only."""
from PySide6.QtWidgets import QWidget, QVBoxLayout, QTextBrowser


GUIDE = """
# RE-FIGS 使用说明

用于 DI-SPA 质谱数据的肽段鉴定与无标记定量。通过左侧导航进入「参数配置」「查看进度」和本说明页面。

## 1　准备输入文件

准备待分析的 **.mzXML** 数据文件和 **.mgf** 参考图谱库，并选择有足够可用空间、可写入的结果文件夹。一次可添加多个数据文件。

## 2　配置分析参数

在「参数配置」页面依次完成以下设置。

### 原始数据导入

- **最终结果存储位置**：选择保存本次分析的文件夹。软件会在其中自动创建批次目录及各文件任务目录。
- **图谱库文件路径**：选择 .mgf 文件。
- **需要分析的 mzXML 文件**：添加待处理文件，开始前移除不需要的条目。
- **CPU 核心数**：设置同时处理的文件数量。数量越大，通常需要更多内存；请根据机器资源选择，实际并行数量不会超过文件数。

### 诱饵库参数

- **使用预计算 decoy 库**：开启后选择包含 library.npy 和 csodiaq_library.npy 的可信目录，例如历史批次的 _shared 文件夹。软件将复用库并跳过建库，无需再选 MGF。策略、距离和 Top 10 沿用预计算时的设置，对应控件禁用；FDR、匹配容差和后续鉴定参数仍可调整。自动寻参仍按开关执行。两个文件必须来自同一次建库，缺失或无法读取时本次分析停止。关闭此选项则恢复从 MGF 建库。
- **策略**：trap-decoy 同时生成诱饵库和反库；decoy 生成反库。
- **entrapment_distance**：诱饵库质量偏移，单位 ppm，通常大于匹配质量容差。
- **decoy_distance**：反库质量偏移，单位 ppm，通常大于 entrapment_distance。
- **FDR 阈值**：设置鉴定筛选阈值，例如 0.01 表示 1%。运行前请确认实际填写值。

### CsoDIAq 参数

- **tolerance**：谱图匹配的质量容差，单位 ppm。
- **is_corrected**：是否启用质量偏移和容差校正。

### RE-FIGS 参数

- **Top 10 的峰**：选择是否仅使用图谱库中强度最高的 10 个峰。
- **start_cycle / end_cycle**：选择扫描循环范围。当前处理包含 start_cycle，不包含 end_cycle；例如 1 到 14 对应第 1–13 个循环。
- **自动搜索 RE-FIGS 参数**：开启后，软件从本批次选取最多 5 个文件搜索参数，再用于正式分析；关闭后使用手动配置的阈值。
- **peaks_shared_limit**：优质训练样本的匹配峰数量阈值。
- **cosine score limit**：余弦相似度阈值。
- **sqrt cosine score**：平方根强度相关相似度阈值。
- **cycle within limit**：优质样本出现的循环数阈值。
- **random seed**：随机种子，用于控制随机步骤；需要比较结果时保持一致。

首次使用可从界面预设参数开始，根据实验数据调整。自动寻参会增加正式分析前的等待时间。

## 3　开始分析并查看进度

配置完成后点击开始执行按钮，进入「查看进度」。软件依次准备图谱库、按需搜索参数，再处理各数据文件。准备和寻参阶段显示忙碌进度；正式处理时显示已结束文件数，失败或取消也计入已结束数量，请同时查看成功与失败数。

文件表显示名称、状态和耗时，可通过上方筛选框按状态查看。超长文件名以中间省略号显示，鼠标悬停查看完整路径。拖动中间分隔条调整文件表与日志区域的宽度。

点击某个文件查看它的实时或历史日志。右侧标题可横向移动光标、选中并复制完整路径。查看历史文件不会被其他文件的新日志切换。

- **搜索**：输入文字后按回车查找，再次回车查找下一处；搜索范围为当前显示日志。
- **复制**：复制当前显示的日志内容。
- **到底部**：回到最新日志；向上滚动查看旧内容时暂停自动跟随。
- **打开日志**：通过系统默认程序查看完整日志。界面最多保留最近 3000 行。
- **打开日志目录 / 打开结果目录**：选中文件后打开该文件对应的任务文件夹；选中准备或寻参阶段时打开共享日志或批次目录。

## 4　终止分析

点击「终止任务」会停止**整个批次**，目前不能只停止单个文件。正在运行的文件标记为「已终止」，未开始的文件标记为「已取消」，已完成文件的结果和日志保留。

终止后软件先等待计算退出，超过 3 秒仍未退出则强制结束。这 3 秒仅在请求终止或关闭窗口后计时，正常计算不会因为耗时长而自动终止。

关闭运行中的窗口也会停止批次并等待清理。中断文件可能保留不完整结果，请以状态和日志为准。任务结束后参数仍保留，调整后可再次运行，新的批次不会覆盖历史目录。

## 5　查找结果和日志

所选结果文件夹中，每次运行生成一个「时间_批次ID」目录。其中每个输入文件对应「序号_文件名」任务文件夹，包含该文件的 CSV 结果和 **task.log**。同名文件通过序号区分。

**_shared** 保存本批次共用的 library.npy、csodiaq_library.npy，以及开启自动寻参时生成的 refigs_params.npy；其 logs 文件夹保存输入参数记录 run.log、建库日志 library.log 和寻参日志 parameter_search.log。

### 常用结果文件

| 文件 | 用途 |
| --- | --- |
| *NoFilter.csv | CsoDIAq 初步谱图匹配结果 |
| *_withFeature_*.csv | 匹配结果及相似度、匹配峰等特征；具体后缀随参数变化 |
| *_pvalue.csv | 评分及统计筛选相关结果 |
| *_LDA_ID.csv | 最终肽段鉴定结果 |
| *_Coeffs.csv | 各扫描中肽段的解卷积系数 |
| *_Quants.csv | 肽段定量结果 |

鉴定结果常见字段包括 peptide（肽段序列）、scan（扫描号）、LDA_score（评分）和 FDR。系数结果中 Coeff 表示肽段在相应扫描中的贡献；定量结果中 Quantify 为定量值，另包含信噪比、峰宽等字段。

## 6　常见操作提示

- **无法开始**：检查是否已选择数据文件、图谱库和可写入的结果目录。
- **文件失败**：选中该文件，打开 task.log 查看完整错误；检查输入格式、参数与可用内存，再调整重跑。
- **准备或寻参失败**：查看对应阶段日志；本次后续文件将被取消。
- **运行慢或内存紧张**：降低 CPU 核心数以减少同时处理的文件数。长时间计算不一定代表程序异常。
- **查找上次结果**：在所选结果目录中打开对应时间的批次文件夹；界面显示当前批次，旧日志仍在磁盘中。
"""


class IntroducationInterface(QWidget):
    def __init__(self, object_name, parent=None):
        super().__init__(parent)
        self.setObjectName(object_name.replace(' ', '-'))
        layout = QVBoxLayout(self)
        layout.setContentsMargins(24, 20, 24, 20)
        self.browser = QTextBrowser(self)
        self.browser.setOpenExternalLinks(True)
        self.browser.document().setDefaultStyleSheet('''
            body { font-family: "Microsoft YaHei UI", "Segoe UI"; font-size: 14px; color: #33465e; }
            h1 { font-size: 27px; color: #1c395b; }
            h2 { font-size: 20px; color: #22659e; margin-top: 26px; }
            h3 { font-size: 16px; color: #243f60; margin-top: 18px; }
            p, li { line-height: 1.65; }
            th { background-color: #edf3fa; }
            td, th { padding: 8px; }
        ''')
        self.browser.setMarkdown(GUIDE)
        self.browser.setStyleSheet('''
            QTextBrowser { background: #ffffff; border: 1px solid #dce4ee; border-radius: 10px; padding: 22px; font-family: "Microsoft YaHei UI", "Segoe UI"; font-size: 14px; color: #33465e; }
            QScrollBar:vertical { background: #f1f4f8; width: 10px; margin: 2px; }
            QScrollBar::handle:vertical { background: #b9c7d7; border-radius: 4px; min-height: 30px; }
            QScrollBar::handle:vertical:hover { background: #7c9bbd; }
            QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0; }
            QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical { background: transparent; }
        ''')
        layout.addWidget(self.browser)
