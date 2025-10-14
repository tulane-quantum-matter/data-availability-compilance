  integer, allocatable :: left(:), right(:)

  integer :: N
  integer :: ml, mr, i, j, k

  N = size(mess);
  allocate(perm(N), left(N/2), right(N - N/2))

  perm = [(i, i=1,N)]
  select case (N)
  case (1)

  case (2)
     if(mess(1) > mess(2)) then
        perm = [2,1]
     end if

  case default
     ml = N/2; mr = N - ml

     left = sortperm(mess(:ml))
     right = sortperm(mess(ml+1:))

     i=1; j=1; k=1
     do while(i<=ml .and. j<=mr)
        if (mess(left(i)) <= mess(right(j)+ml)) then
           perm(k) = left(i)
           i=i+1
        else
           perm(k) = right(j)+ml
           j=j+1
        end if
        k=k+1
     end do
     if(i<=ml) then
        perm(k:) = left(i:)
     else
        perm(k:) = right(j:)+ml
     end if

  end select

  deallocate(left, right)
