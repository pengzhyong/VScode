#include "predict.h"
#include "io.h"
#include "stdio.h"

int main()
{
    print_time("Begin");
	char *data[MAX_DATA_NUM];
    char *info[MAX_INFO_NUM];
	int data_line_num;
    int info_line_num;


	char *data_file = "data/TrainData_2015.1.1_2015.2.19.txt";

    data_line_num = read_file(data, MAX_DATA_NUM, data_file);

    printf("data file line num is :%d \n", data_line_num);
    if (data_line_num == 0)
    {
        printf("Please input valid data file.\n");
        return -1;
    }
	
	char *input_file = "data/input_5flavors_cpu_7days(1).txt";

    info_line_num = read_file(info, MAX_INFO_NUM, input_file);

    printf("input file line num is :%d \n", info_line_num);
    if (info_line_num == 0)
    {
        printf("Please input valid info file.\n");
        return -1;
    }

	char *output_file = "out_file.txt";

    predict_server(info, data, data_line_num, output_file);

    release_buff(info, info_line_num);
	release_buff(data, data_line_num);

    print_time("End");

	return 0;
}

