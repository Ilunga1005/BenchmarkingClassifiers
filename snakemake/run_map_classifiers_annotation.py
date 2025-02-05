import argparse
import logging # 用于记录日志。
import sys # 用于与系统交互，如退出程序。

import yaml # 用于处理YAML格式的配置文件。
import snakemake #执行SnakeMake工作流。
from pathlib import Path
from rich.logging import RichHandler # 用于更美观的日志输出（通过 rich 库）。
from typing import Dict, Any, Optional, Sequence
from datetime import datetime


def run_snakemake(dir_working: Path, snakefile: Path, dry_run: bool, config_data: Dict[str, Any]) -> None:
    """
    Runs snakemake.
    :param dir_working: Directory in which the *.smk will be executed
    :param snakefile: Snakefile
    :param dry_run: If True, will execute the SnakeMake command with parameter -n
    :param config_data: Snakemake configuration data
    :return:None
    """
    # Create working directory
    if not dir_working.exists():
        dir_working.mkdir(parents=True)

    # Create snakemake config
    config_path = (dir_working / 'snakemake_config.yml').resolve()
    with config_path.open('w') as handle:
        yaml.dump(config_data, handle)
    logger.info(f"Snakemake config created: {config_path}")

    snakemake.snakemake(snakefile=snakefile,
                        configfiles=[config_path],
                        workdir=dir_working,
                        stats=dir_working / ("stats_" + datetime.now().strftime("%Y_%m_%d_%H_%M_%S") + ".json"),
                        printreason=True,
                        printshellcmds=True,
                        cores=32,
                        dryrun=dry_run
                        )


class Classifying(object):
    """
    This class contains is used to perform the classifying .
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None: # 类的构造函数，它在创建类的实例时自动调用，用于初始化对象的属性。 
        """
        Initializes the classifying context instance.
        args: 为一个参数，表示传递给类的参数列表。
        Optional[Sequence[str]] 是类型注解，表示 args 可以是 Sequence[str] 类型（例如列表或元组），也可以是 None。
        = None 是默认值，如果调用时没有传递 args，则默认为 None。
        -> None：这是返回类型注解，表示该方法没有返回值（返回 None）。
        """
        self._args = Classifying._parse_arguments(args) # 调用 _parse_arguments 方法解析命令行参数，并将结果赋值给 self._args。
        """
        _parse_arguments(args), 单下划线 _ 开头的属性通常表示它是受保护的（protected），即不建议外部代码直接访问或修改。
        """

    @staticmethod
    # @staticmethod 是一个装饰器，用于定义一个静态方法。静态方法不需要访问实例属性或方法，因此不需要 self 参数。
    #类似于一个普通函数，但是定义在了类的命名空间中
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Arguments
        :return: Parsed arguments
        返回类型是 argparse.Namespace，表示解析后的命令行参数。
        """
        parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('--input',
                            help="FASTA/Q input file",
                            required=True,
                            type=lambda x: Path(x).resolve() # 将输入的路径字符串转换为 Path 对象，并解析为绝对路径。
                            )
        parser.add_argument('--truth', ## ？？？？
                            help='Location to the ground truth file.',
                            required=True
                            )

        parser.add_argument('--output',
                            help='Output of snakemake',
                            required=True,
                            type=lambda x: Path(x).resolve())
        parser.add_argument('--config_file',
                            help='YAML file with location of tools and database.\n'
                                 'Template can be found under "/snakemake/config/classifiers.yml.template"',
                            required=True
                            )
        parser.add_argument('--classifier',
                            nargs='+',
                            default=['bracken', 'ccmetagen', 'centrifuge', 'kaiju', 'kma', 'kraken2', 'metaphlan',
                                     # 'mmseqs2',
                                     'motus'],
                            help='Which classifier(s) to use. \n'
                                 'Default choices are %(default)s; space separated.\n'
                                 'If not set, all default choices will be used.',
                            metavar='' # metavar=''：在帮助信息中隐藏参数的元变量名称。
                            )
        parser.add_argument('--dir_working',
                            help='Working directory. \n'
                                 'If not set, it will be {output path} + "snakemake_conf".',
                            default=None
                            )
        parser.add_argument('--dry_run', help="Run the SnakeMake file in Dry Run mode", action='store_true')
        parser.add_argument('--no_filter', help='Disable filtering of reads (length > 1000 and quality > 7)', action='store_true')

        parser_args = parser.parse_args(args)

        return parser_args

    def _check_input_arguments(self):
        """
        Check if input file really exists and arguments make sense
        """
        if not self._args.input.exists():
            logger.error(f'Input file {self._args.input} does not exist!')
            sys.exit(1)

        mutual_inclusive = {'ccmetagen': 'kma', 'bracken': 'kraken2'} # 定义了一个字典，其中键是分类器的名称，值是必须同时使用的分类器的名称。定义了分类器之间的依赖关系。
        for name, dependency in mutual_inclusive.items():
            if name in self._args.classifier and dependency not in self._args.classifier:
                logger.error(f"Error: '{dependency}' is required when '{name}' is selected.")
                sys.exit(1)

        # Set working directory to output directory if not specified
        if self._args.dir_working is None:
            self._args.dir_working = self._args.output / 'snakemake_conf'
            logger.info(f"Working DIR set to {self._args.dir_working}")

    def __get_config_data(self) -> Dict[str, Any]:
        """
        Get the snakemake configuration data.
        :return: Configuration data
        """

        with open(self._args.config_file, 'r') as f:
            locations = yaml.load(f, Loader=yaml.FullLoader)

        return {
            'output': str(self._args.output),
            'input': {'fastq': str(self._args.input)},
            'no_filter': self._args.no_filter,
            'locations': locations,
            'classifiers': self._args.classifier,
            'truth_file': self._args.truth,
            'python_scripts': str(Path(__file__).parent.resolve() / 'scripts')
        }

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        # Check if input exists
        self._check_input_arguments()

        run_snakemake(
            self._args.dir_working,
            (Path(__file__).parent / 'smk' / 'map_classifiers.smk').resolve(),
            self._args.dry_run,
            self.__get_config_data())


if __name__ == '__main__': # 确保只有当脚本作为主程序执行时，以下代码才会被执行。

    logger = logging.getLogger(__name__) # 使用__name__的名字作为日志记录器的名字
    logger.setLevel(logging.DEBUG) # 设置日志记录器的日志级别为 DEBUG，意味着所有 DEBUG 级别及以上的日志信息（例如：INFO、WARNING、ERROR 和 CRITICAL）都会被记录和输出。
    if not logger.handlers:
        # create console handler and set level to debug
        ch = logging.StreamHandler() # 创建一个控制台处理器，用于将日志信息输出到控制台。即将日志信息输出到终端
        ch.setLevel(logging.DEBUG) # 设置日志级别

        # create formatter
        formatter = logging.Formatter('%(levelname)s: %(message)s') # 日志输出格式，包括日志级别和日志信息。
        # add formatter to ch
        ch.setFormatter(formatter) # 将格式化器添加到控制台处理器中，这样控制台处理器在输出日志信息时就会使用这个格式化器。
        # add ch to logger
        # logger.addHandler(ch)
        logger.addHandler(RichHandler(show_path=False, level=1)) # 将 RichHandler 作为日志处理器添加到日志记录器中，这样日志消息将会通过这个处理器输出

    classify = Classifying() # 创建一个 Classifying 类的实例
    classify.run()
