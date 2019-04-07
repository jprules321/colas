import {Component, OnInit} from '@angular/core';
import {ElectronService} from '../../providers/electron.service';
import {HttpClient} from '@angular/common/http';
import {Observable} from 'rxjs';

@Component({
    selector: 'app-home',
    templateUrl: './home.component.html',
    styleUrls: ['./home.component.scss']
})
export class HomeComponent implements OnInit {

    options;
    datas: any[] = [];


    constructor(private electron: ElectronService, private http: HttpClient) {
    }

    public getGribDatas = (event) => {
        const exe = 'src\\wgrib2\\wgrib2.exe';
        const file = `${event.target.files[0].path}`;
        const options = this.options ? this.options : '';
        for (let i = 0; i < 500;){
            setTimeout(() => {
                this.electron.cmd.get(
                    `${exe} ${file} ${options}`,
                    function(err, data, stderr){
                        console.log(data || stderr || err);
                        i++;
                    }
                );
            }, 10)
        }
    };



    ngOnInit() {
    }

}
