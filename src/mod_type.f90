module mod_type
  use iso_fortran_env, only: dp => real64,&
    li => int64,stdout => output_unit,&
    si => int32
  implicit none
  real(dp), parameter :: SQR2 = 1.4142135623730950488d0     !< 根号2

  private
  public :: Frame, Data

  type :: Frame
    character(:),allocatable :: prefix   !< 每一帧文件的前缀
    integer(si) :: frames                !< 帧数
    character(:),allocatable :: suffix   !< 文件的后缀
    character(:),allocatable :: filename !< 上面三个合起来就是完整的文件名称
    character(:),allocatable :: property !< 获取dump文件中输出了哪些类型的property
    integer(si) :: num_of_atoms          !< 获得总共有多少个原子
    integer(si),allocatable :: ids(:)    !< 每个原子的id数组
    integer(si),allocatable :: types(:)  !< 每个原子的类型
    real(dp),allocatable :: datas(:,:)   !< 每个原子的具体物理属性数值
    real(dp) :: xu,xl,yu,yl,zu,zl        !< l表示某个方向的下限,u表示某个方向的上限
  end type Frame

  type Data
    integer(si) :: atom_types                             !< 多少种原子
    integer(si) :: natoms                                 !< 一共多少个原子
    integer(si), allocatable, dimension(:) :: type_array  !< 类型数组
    real(dp), allocatable, dimension(:, :) :: coord       !< 记录原子的坐标和类型,第一列为类型, 二到四列为坐标
    real(dp), dimension(3, 2) :: box                      !< 记录盒子的极限坐标
  contains
    procedure, public, pass(self) :: WCSRO
  end type Data

  interface Data
    module procedure :: Data_constructor                  !< 覆盖掉Data这个类的默认构造器
  end interface Data

  interface Frame
    module procedure :: Frame_constructor
  end interface Frame

contains

  subroutine WCSRO(self, ltype, lattice)
    !! 计算模型的WCSRO系数, 默认坐标按照xyz排序
    class(Data), intent(in) :: self         !< 包含了原子位置和盒子坐标的数据集合
    character(len=*), intent(in) :: ltype   !< 晶格类型
    real(dp), intent(in) :: lattice         !< 模型的晶格常数
    real(dp) :: rCut                        !< 截断半径
    real(dp) :: lx, ly, lz                  !< 盒子半径
    real(dp) :: dR                          !< 原子之间的距离
    real(dp), allocatable :: ratio(:)       !< 存储每个原子的百分比
    integer, allocatable :: counts(:)       !< 存储每种原子的近邻原子数量
    real(dp), allocatable, dimension(:, :) :: WCP  !< 用于存放用来计算wallon-colley parameter的数组
    integer :: i, j                                !< 循环变量

    !! 根据传入的晶格类型和晶格常数计算截断半径
    select case(trim(ltype))
     case ('fcc')
      rCut = 0.25d0 * (SQR2 + 2.d0) * lattice
     case('bcc')
      rCut = 0.5d0 * (SQR2 + 1.d0) * lattice
     case default
      stop "无法识别的晶格类型, 本程序仅支持fcc和bcc两种晶格类型的计算, 请按正确格式输入."
    end select

    rCut = rCut * rCut    !! 这里我们将对比的是原子之间距离的平方和截断半径距离的平方, 以避免不必要的开方操作
    !! 计算盒子的几个边长
    lx = self%box(1, 2) - self%box(1, 1)
    ly = self%box(2, 2) - self%box(2, 1)
    lz = self%box(3, 2) - self%box(3, 1)

    !! 计算每个原子的百分比
    allocate(ratio(self%atom_types), counts(self%atom_types))
    ratio = 0.d0; counts = 0
    do i = 1, self%atom_types
      ratio(i) = 1.d0 - count(self%type_array /= i) / real(self%natoms, kind=dp)
    end do

    !! 分配内存,并初始化
    allocate(WCP(self%atom_types, self%atom_types))
    WCP = 0.d0

    !! 循环所有原子对, 填补WCP数组, 调用openMP并行
    !$OMP parallel do default(firstprivate) shared(self, WCP, rCut, lx, ly, lz, counts) schedule(dynamic, 2)
    do i = 1, self%natoms - 1
      do j = i + 1, self%natoms
        dR = Distance(self, i, j, lx, ly, lz)   !! 计算2个原子之间的距离, 考虑周期性边界条件, 在三个方向上
        if (dR < rCut) then
          !! 注意如果是有条件的对同一个原子进行数学操作, 一定要加入!$OMP critical保证不出错
          !$OMP critical
          counts(self%type_array(i)) = counts(self%type_array(i)) + 1
          counts(self%type_array(j)) = counts(self%type_array(j)) + 1
          WCP(self%type_array(i), self%type_array(j)) = WCP(self%type_array(i), self%type_array(j)) + 1.D0
          WCP(self%type_array(j), self%type_array(i)) = WCP(self%type_array(j), self%type_array(i)) + 1.D0
          !$OMP end critical
        end if
      end do
    end do
    !$OMP end parallel do
    !! 在进行了最耗时的计算之后, 通过循环计算更新WCP
    do i = 1, self%atom_types
      do j = i, self%atom_types
        WCP(i, j) = 1.d0 - WCP(i, j) / (ratio(j) * counts(i))   !! 根据定义表达式计算WCP参数
        if (i /= j) WCP(j, i) =  WCP(i, j)                      !! 这里其实很容易从数学表达式推出WCP_mn = WCP_nm
      end do
    end do
    !! 格式化输出WCP数组看看
    call fprint_WCP(WCP, self%atom_types)

  end subroutine WCSRO

  subroutine fprint_WCP(WCP, Num)
    !! 将计算出来的WCP格式化输出
    real(dp), intent(in), dimension(:, :) :: WCP
    integer, intent(in) :: Num
    integer :: i


    do i = Num, 1, -1
      write(stdout, '(I5, 3X, *(F5.2, 3X))') i, WCP(i, :)
    end do
    write(stdout, '(T9, *(I5, 3X))') (i, i = 1, Num)

  end subroutine

  pure function Distance(Atom, i, j, lx, ly, lz) result(res)
    type(Data), intent(in) :: Atom
    real(dp), intent(in) :: lx, ly, lz
    integer, intent(in) :: i, j
    real(dp) :: res, dx, dy, dz

    dx = Atom%coord(i, 1) - Atom%coord(j, 1)
    dx = dx - nint(dx / lx) * lx
    dy = Atom%coord(i, 2) - Atom%coord(j, 2)
    dy = dy - nint(dy / ly) * ly
    dz = Atom%coord(i, 3) - Atom%coord(j, 3)
    dz = dz - nint(dz / lz) * lz

    res = dx * dx + dy * dy + dz * dz

  end function Distance

  function Data_constructor(filename) result(res)
    character(len=*), intent(in) :: filename  !< 文��名称
    type(Data) :: res
    integer :: fileid, io_flag, i, j
    logical :: is_exist

    inquire(file=filename, exist=is_exist)
    if (.not. is_exist) then
      write(stdout, '(A, A, A)') 'file ', filename, ' does not exist!'
      stop 
    end if

    open(newunit=fileid, file=filename, action='read')
    do i = 1, 2
      read(fileid, *, iostat=io_flag)
      if (io_flag /= 0) then
        write(stdout, '(A, A)') 'Error occur when reading datafile ', filename
        stop
      end if
    end do
    read(fileid, *, iostat=io_flag) res%natoms
    if (io_flag /= 0) then
      write(stdout, '(A, A)') 'Error occur when reading numbers of atoms '
      stop
    end if
    read(fileid, *, iostat=io_flag)
    read(fileid, *, iostat=io_flag) res%atom_types

    if (io_flag /= 0) then
      write(stdout, '(A, A)') 'Error occur when reading numbers of atom types '
      stop
    end if

    read(fileid, *, iostat=io_flag)     !! 所有的空read都是在跳过注释行, 一个read跳过一行

    do i = 1, 3
      read(fileid, *, iostat=io_flag) res%box(i, 1), res%box(i, 2)
      if (io_flag /= 0) then
        write(stdout, '(A, A)') 'Error occur when reading numbers of coordination of box '
        stop
      end if
    end do

    do i = 1, 3
      read(fileid, *, iostat=io_flag)
    end do

    allocate(res%type_array(res%natoms), res%coord(res%natoms, 3))

    do i = 1, res%natoms
      read(fileid, *, iostat=io_flag) j, res%type_array(i), res%coord(i, :)
      if (io_flag /= 0) stop 'Error occur when reading information of coordination.'
    end do
    close(fileid)

  end function Data_constructor

  function num_records(filename) result(res)
    !! 此函���用于读取文件里面有多少行数据
    character(len=*),intent(in) :: filename
    integer(li) :: res,temp
    logical :: alive

    res = 0   ! 初始为0

    ! 查询文件是否存在
    inquire(file=filename,exist=alive)

    ! 如果文件不存在就终止程序运行
    if (.not. alive) then
      write(stdout,*) trim(filename),' does not exist !'
      stop
    end if
    open(newunit=temp,file=filename,action='read')

    ! 反��读取数据,直������������������������一行,���至1号进行������步处理
    do
      read(temp,*,end=1)
      res = res+1
    end do

1   close(temp)   ! 关���文件
    return

  end function num_records

  type(frame) function Frame_Constructor(prefix,frames,suffix) result(self)
    character(*),intent(in) :: prefix
    integer(si), intent(in) :: frames
    character(*),intent(in) :: suffix
    character(20) :: num_frame
    character(:),allocatable :: property_text
    integer(si) :: num_atoms
    integer(si) :: property
    integer(si),allocatable :: ids(:)
    integer(si),allocatable :: types(:)
    real(dp),allocatable :: datas(:,:)
    real(dp) :: xu,xl,yu,yl,zu,zl

    ! 赋予初值
    self%prefix = prefix
    self%frames = frames
    self%suffix = suffix

    ! 首先将整数frames改成字符串变量
    write(num_frame,"(G0)") frames

    ! 获得完整的文件名称
    self%filename = trim(prefix)//trim(num_frame)//trim(suffix)

    ! 读取文件,������得关键信������
    call read_dump(self%filename,num_atoms,&
      datas,ids,types,property,property_text,&
      xu,xl,yu,yl,zu,zl)

    ! 利用读取到的文件对frame类������面的数组分配内����
    allocate(self%ids(num_atoms),self%types(num_atoms),&
      self%datas(num_atoms,property))

    ! ���配内存之后将读取到的数据��值给frame��中的数据域
    self%num_of_atoms = num_atoms
    self%ids = ids
    self%types = types
    self%datas = datas
    self%property = property_text
    self%xl = xl;   self%xu = xu
    self%yl = yl;   self%yu = yu
    self%zl = yl;   self%zu = yu

  end function Frame_Constructor

  subroutine read_dump(filename,num,dataline,ids,&
    types,property,property_text,&
    xu,xl,yu,yl,zu,zl)
    !! 对于intent类型为out的动态数组,传��的时候��自动deallocate
    character(len=*),intent(in) :: filename
    character(:),allocatable,intent(out) :: property_text
    integer(si),intent(out) :: num
    real(dp),allocatable,intent(out) :: dataline(:,:)
    integer(si),allocatable,intent(out) :: ids(:)
    integer(si),allocatable,intent(out) :: types(:)
    integer(si),intent(out) :: property
    integer(si) :: fileunit,i
    character(len=4096) :: textline
    logical :: alive
    real(dp),intent(out) :: xu,xl,yu,yl,zu,zl

    ! 查询文件是否存��
    inquire(file=filename,exist=alive)

    ! 如果文件�����在就终止程序运行
    if (.not. alive) then
      write(stdout,*) trim(filename),' does not exist !'
      stop
    end if

    open(newunit=fileunit,file=trim(filename),&
      action='read',status='old')

    ! 跳��������行文本��
    do i = 1,3
      read(fileunit,*)
    end do

    ! 从第四行的数据��获取原子总数
    read(fileunit,*) num
    read(fileunit,*)              ! 跳过这一行文本
    ! 分配内存
    allocate(ids(num),types(num))

    ! 跳�����前4��文本行
    read(fileunit,*) xl,xu
    read(fileunit,*) yl,yu
    read(fileunit,*) zl,zu


    ! 从第九行的数据中获得输�����具体物理属性,并判断是否符合格式规范
    read(fileunit,'(A)') textline
    if (textline(13:20) /= 'id type') then
      write(stdout,"(A,A)") 'Error occuring while reading: ',trim(filename)
      write(stdout,"(A)") 'The 1st and 2nd Column must be [id type] !'
      write(stdout,'(A)') 'Please correct the format of datafile.'
      stop
    end if
    property_text = trim(textline)
    property = read_property(textline)

    ! 分配����存
    allocate(dataline(num,property))

    do i = 1, num
      read(fileunit,*) ids(i),types(i),dataline(i,:)
    end do

    ! 关闭文件
    close(fileunit)

  end subroutine read_dump

  pure function read_property(text) result (res)
    character(len=*),intent(in) :: text
    integer(si) :: res,i,temp

    ! temp为text字符串的��际长度
    temp = len_trim(text); res = 0

    ! 遍���整���text字符串,如果��现空格,res就���一
    do i = 1, temp
      if (text(i:i) == ' ') res = res + 1
    end do

    res = res - 3  ! 减二之后才是真正的输��的property���数

  end function read_property

end module mod_type
