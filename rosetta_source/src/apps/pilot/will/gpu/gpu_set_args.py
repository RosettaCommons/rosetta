t = ("cl_mem","   int","  uint"," float","")

done = set()
for i in range(len(t)-1):
    for j in range(len(t)):
        for k in range(len(t)):
            for l in range(len(t)):
                #for m in range(5):
                    args = ", "+t[i]+" a1"
                    if t[j]: 
                        args += ", "+t[j]+" a2"
                        if t[k]:
                            args += ", "+t[k]+" a3"
                            if t[l]:
                                args += ", "+t[l]+" a4"
                                #if t[m]:
                                #    args += ", "+t[m]+" a5"
                    if args in done: continue
                    print "  void setargs(std::string kname"+args+") {\n    std::string const s(\"Error: Failed to set kernel argument \"); err_ = 0;"
                    if 1:                      print "    err_=clSetKernelArg(kernels_[kname],0,sizeof("+t[i]+"),&a1); if(err_!=CL_SUCCESS) handle_error(s+\"1 of "+str(args.count(","))+" type "+t[i]+" for \"+kname+\"! ERR=\"+str(err_));"
                    if t[j]:                   print "    err_=clSetKernelArg(kernels_[kname],1,sizeof("+t[j]+"),&a2); if(err_!=CL_SUCCESS) handle_error(s+\"2 of "+str(args.count(","))+" type "+t[j]+" for \"+kname+\"! ERR=\"+str(err_));"
                    if t[j] and t[k]:          print "    err_=clSetKernelArg(kernels_[kname],2,sizeof("+t[k]+"),&a3); if(err_!=CL_SUCCESS) handle_error(s+\"3 of "+str(args.count(","))+" type "+t[k]+" for \"+kname+\"! ERR=\"+str(err_));"
                    if t[j] and t[k] and t[l]: print "    err_=clSetKernelArg(kernels_[kname],3,sizeof("+t[l]+"),&a4); if(err_!=CL_SUCCESS) handle_error(s+\"4 of "+str(args.count(","))+" type "+t[l]+" for \"+kname+\"! ERR=\"+str(err_));"
                    #print "    err_=clSetKernelArg(kernels_[kname],4,sizeof("+t[m]+"),&a5); if(err_!=CL_SUCCESS) handle_error(s+\"5 of "+str(args.count(","))+" type "+t[m]+" for \"+kname+\"! ERR=\"+str(err_));"
                    print "  }"
                    done.add(args)
