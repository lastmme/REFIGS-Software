# RE-FIGS-Software

## 源代码-环境配置

请按照以下步骤配置环境：

1. **使用 `environment.yml` 文件创建 Conda 环境**：

    ```bash
    conda env create -f environment.yml
    ```

    这将根据 `environment.yml` 文件中的配置创建一个新的 Conda 环境。

2. **使用 `requirements.txt` 文件安装依赖**：

    激活你的 Conda 环境：

    ```bash
    conda activate refigs
    ```

    然后使用 `pip` 安装依赖：

    ```bash
    pip install -r requirements.txt
    ```

    这将根据 `requirements.txt` 文件中的配置安装所有必要的 Python 包。

    3. **注意**：

        请注意，`environment.yml` 文件已经包含了 `requirements.txt` 文件中的配置。因此，通常情况下，你只需要使用 `environment.yml` 文件来创建和配置环境，而不需要单独运行 `pip install -r requirements.txt`。

## 软件安装地址

[RE-FIGS-Software](https://github.com/lastmme/REFIGS-Software/releases/tag/test)
