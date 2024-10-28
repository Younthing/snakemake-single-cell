# 规则：hello_world
# 功能：
#   生成一个简单的文本文件，内容为"Hello World!"消息
rule hello_world:
    output:
        # 输出文件路径，统一命名为hello_world.txt
        "results/example/hello_world.txt",
    shell:
        # 使用shell命令创建文件并写入消息
        "echo 'Hello World!' > {output}"


# 规则：hello_with_params
# 功能：
#   生成一个文本文件，包含重复的自定义消息
# 参数：
#   message：要写入的自定义消息
#   repeat_count：消息重复的次数
rule hello_with_params:
    output:
        # 输出文件，统一命名为hello_with_params.txt
        "results/example/hello_with_params.txt",
    params:
        # 自定义参数
        message="Hello from params!",
        repeat_count=3,
    shell:
        # 使用shell脚本将消息重复写入文件
        """
        for i in $(seq {params.repeat_count}); do
            echo "{params.message}" >> {output}
        done
        """


# 规则：hello_with_input
# 功能：
#   处理一个输入的文本文件，并在其后附加时间戳
# 输入：
#   txt：要处理的输入文件，来自hello_world的输出
# 输出：
#   txt：处理后的文本文件
rule hello_with_input:
    input:
        # 指定输入文件路径，与hello_world的输出一致
        "results/example/hello_world.txt",
    output:
        # 输出文件路径，统一命名为hello_with_input.txt
        "results/example/hello_with_input.txt",
    shell:
        # 读取输入文件内容，复制到输出文件，并添加时间戳
        """
        cat {input} > {output}
        echo "Processed at: $(date)" >> {output}
        """


# 规则：hello_with_config
# 功能：
#   使用配置文件中的参数生成文本文件
# 配置参数：
#   输出目录：config["paths"]["output_dir"]
#   消息内容：config["parameters"]["message"]
#   线程数：config["parameters"]["threads"]

import os
from pathlib import Path


# 加载配置文件
configfile: "config/config_example.yaml"


rule hello_with_config:
    output:
        # 使用配置文件中的路径拼接输出文件路径，统一命名为hello_with_config.txt
        os.path.join(config["paths"]["output_dir"], "hello_with_config.txt"),
    params:
        # 从配置文件中获取自定义消息
        message=config["parameters"]["message"],
    threads: config["parameters"]["threads"]  # 使用配置文件中的线程数
    shell:
        # 将消息写入输出文件
        "echo '{params.message}' > {output}"


# 规则：example_base
# 功能：
#   汇总所有前面规则生成的输出文件
# 输入：
#   所有之前规则生成的输出文件
rule example_base:
    input:
        # 列出所有需要汇总的文件路径
        "results/example/hello_world.txt",
        "results/example/hello_with_params.txt",
        "results/example/hello_with_input.txt",
        os.path.join(config["paths"]["output_dir"], "hello_with_config.txt"),
