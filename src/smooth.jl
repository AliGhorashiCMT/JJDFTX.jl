function smooth(y; win_len=40, win_method=2)
if win_len%2 == 0
    win_len+=1 # only use odd numbers
  end
  if win_method == 1
    w=ones(win_len)
  elseif win_method==2
    w=DSP.hanning(win_len)
elseif win_method==3
w=DSP.hamming(win_len)
end
if win_len<3
return y
elseif length(y)<win_len
return y
else
y_new = [2*y[1] .- reverse(y[1:win_len],1); y[:]; 2*y[end] .- reverse(y[end-win_len:end],1)]
y_smooth = conv(y_new, w/sum(w))
ind = floor(Int, 1.5*win_len)
return y_smooth[1+ind:end-ind-1]
end
    
end 
