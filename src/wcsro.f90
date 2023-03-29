program main
  use iso_fortran_env, only: dp => real64, stdout => output_unit
  use mod_type, only: Data
  implicit none
  character(len=80) :: filename     !< 模型文件的名称
  character(len=80) :: ltype        !< 晶格类型, 目前只有fcc和bcc两个选项
  character(len=80) :: inpurchar    !< 用户输入的第一个字符串
  real(dp) :: lattice               !< 晶格常数
  logical :: is_exist               !< 判断文件是否存在的逻辑变量
  integer :: io_stat                !< 读取过程中, 如果io_stat不为0, 说明出了一些问题
  type(Data) :: Model               !< 存放模型文件中的信息的集合
  namelist /inputparameter/ filename, ltype, lattice    !< 将上述变量放到namelist的集合中
  integer :: i, fileid

  call get_command_argument(NUMBER=1, value=inpurchar)  !! 从命令行中读取用户的第一个输入

  if (len_trim(inpurchar) == 0) stop "Please input filename contains namelist!"
  inquire(file=trim(inpurchar), exist=is_exist)         !! 查询文件是否存在

  if (.not. is_exist) then                              !! 如果文件不存在就报错并退出
    write(stdout, "(A, A, A)") 'File ', trim(inpurchar), " doesn't exist!"
    stop
  end if 

  open(newunit=fileid, file=trim(inpurchar), action='read')   !! 打开文件
  read(UNIT=fileid, NML=inputparameter, iostat=io_stat)       !! 以namelist的标准读取数据

  !! 此时我们已经获得了计算WCP的一切信息, 将其投入子程序中计算即可
  Model = Data(trim(filename))                  !! 初始化数据集合
  call Model%WCSRO(trim(ltype), lattice)        !! 计算WCP参数

end program main
