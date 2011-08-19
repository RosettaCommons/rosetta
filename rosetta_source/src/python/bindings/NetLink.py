#!/usr/bin/env python
# :noTabs=true:

import  socket, time #, zlib
from array import array


class PR_UDPServer:
    def __init__(self, udp_port=65000, udp_ip = '127.0.0.1'):
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        self.socket.bind( (udp_ip, udp_port) )
        self.buf = {}
        self.last_cleanup_time = time.time()



    def listen(self):
        data, addr = self.socket.recvfrom(1024*64)  # 64k buffer
        #print 'Got Message from %s' % str(addr), len(data)

        packet_id = data[:16+2]
        counts = array('H', data[18:22])  # should be 2 short integer

        #print 'Packet count info:', counts
        if counts[1] == 1:  # only one messgage in pack...
            return array('c', data[22:]) #bz2.decompress(data[22:])
        else:
            if packet_id not in self.buf: self.buf[packet_id] = [0., {}]

            c, d = self.buf[packet_id]
            d[counts[0]] = data[22:]
            self.buf[packet_id][0] = time.time()

            # now, lets check if we can find all the pieces of this message...
            if len(d) == counts[1]:  # yes, they all here...
                #print 'Asseblinbg messge from %s pieces...' % counts[1]
                m = array('c')
                for i in range(counts[1]):
                    m.extend( d[i] )
                del self.buf[packet_id]

                #print 'Messge is:', len(m), m, d
                #print 'leftover buffer len:', len(self.buf)
                return m #bz2.decompress(m)

            else:
                # there is no ready-to-return packets... however lets check if buffer can be cleaned up...
                # anything older then 10second should be discarded...
                current_time = time.time()
                if current_time - self.last_cleanup_time > 2.: # cleaning up every 2s

                    for k in self.buf.keys():
                        if current_time - self.buf[k][0] > 10.0 : del self.buf[k]

                return None




if __name__ == "__main__":
    from cStringIO import StringIO
    import gzip


    ps = PR_UDPServer()
    while True:
        r = ps.listen()
        if r:
            #r = zlib.decompress(r)
            rn = gzip.GzipFile('', 'r', 0, StringIO(r)).read()
            print rn, len(rn), len(r)




if __name__ == "__main__":
    import time

    while True:

        #message = 'pdb ' + '\x05decoy' + '.' * 1024*64  + ' Mda...'*10245 + 'r'*1024*1024*1
        message = 'pdb ' + '\x05decoy' + file('test_in.pdb').read()
        us = PySocketClient()
        us.sendMessage(message)

        #time.sleep(1)


