import unittest

from run_spaceranger import build_spaceranger_command


class SpaceRangerCommandTests(unittest.TestCase):
    # 功能：验证原始 v4 命令中的关键参数能由 YAML 模板生成。
    def test_build_command(self):
        command = build_spaceranger_command(
            {
                "id": "sample1",
                "transcriptome": "/ref",
                "fastqs": ["/fq1", "/fq2"],
                "image": "/image.tif",
                "create_bam": False,
            }
        )
        self.assertIn("--id=sample1", command)
        self.assertIn("--fastqs=/fq1,/fq2", command)
        self.assertIn("--create-bam=false", command)


if __name__ == "__main__":
    unittest.main()
