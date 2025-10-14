  integer :: N, ml, mr, i, j, k
  N = size(mess)
  allocate(left(N/2), right(N - N/2))

  select case (N)
  case (1)

  case (2)
     if(mess(1) > mess(2)) then
        mess(1:2) = mess(2:1:-1)
     end if

  case default
     ml = N/2; mr = N - ml
     left = mess(:ml)
     right = mess(ml+1:)
     call sort(left)
     call sort(right)

     i=1; j=1; k=1
     do while(i<=ml .and. j<=mr)
        if(left(i) <= right(j)) then
           mess(k) = left(i)
           i=i+1
        else
           mess(k) = right(j)
           j=j+1
        end if
        k=k+1
     end do
     if(i<=ml) then
        mess(k:) = left(i:)
     else
        mess(k:) = right(j:)
     end if
  end select

  deallocate(left,right)
