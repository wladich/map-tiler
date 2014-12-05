# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
MIT License
===========
Copyright (c) 2013 Sergej Orlov
 
Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:
 
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE."""

import multiprocessing
import traceback
import sys
import signal
import itertools

class ChildException(Exception):
    pass
    
def mpimap_wrapper((func, args, kwargs)):
    result = {'error': None}
    try:
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        result['value'] = func(*args, **kwargs)
    except Exception:
        err_cls, err, tb = sys.exc_info()
        result['error'] = err
        result['traceback'] = ''.join(traceback.format_tb(tb))
    return result
    
def mpstarimap(func, job, **kwargs):
    no_mp = kwargs.pop('__nomp', False)
    if no_mp:
        imap = itertools.imap
    else:
        pool = multiprocessing.Pool()
        imap = pool.imap
    job = ((func, args, kwargs) for args in job)
    for result in imap(mpimap_wrapper, job):
        error = result['error']
        if error:
            raise ChildException('%r\n%s' % (error, result['traceback']))
        yield result['value']
        
def mpimap(func, job, **kwargs):
    job = ((x,) for x in job )
    return mpstarimap(func, job, **kwargs)

if __name__ == '__main__':
    import time

    def iterate_jobs():
        for i in xrange(100):
            yield i, i*i
    
    def process(x, y, z):
        time.sleep(0.1)
        if x > 20:
            raise ValueError(x)
        return x + y + z
    
    n3 = 22

    for result in mpimap(process, iterate_jobs(), z=n3):
        print result