from typing import List
import numpy as np
from scipy import sparse
import osqp
import astra

class Identified:
    id_: int
    name_: str
    index_: int

    def __init__(self, index: int, id: int, name: str):
        self.index_ = index
        self.id_ = id
        self.name_ = name

    def verbal_name(self):
        return f"{self.id_:8} [{self.name_[:20]:20}]"

class Load(Identified):
    pbase_: float
    pmin_: float
    pmax_: float
    node_id_: int
    const: bool

    def __init__(
        self,
        index: int,
        id: int,
        name: str,
        node_id: int,
        pbase: float,
        pmin: float,
        pmax: float,
        const: bool,
    ):
        super().__init__(index, id, name)
        self.node_id_ = node_id
        self.pbase_ = pbase
        self.pmin_ = pmin
        self.pmax_ = pmax
        self.const_ = const


class LoadMap:
    loads_: List[Load]
    load_keys_: dict[int, int]

    def __init__(self):
        self.loads_ = []
        self.load_keys_ = {}

    def check_assert(self):
        for load in self.loads_:
            if load.const_ or load.pmax_ - load.pmin_ < 1e-6:
                raise Exception(f"Нагрузка {load.verbal_name()}")


class Territory(Identified):
    pbase_: float
    pbase_by_loads_: float
    uncontrolled_load_: float
    lrc_load_: float
    node_indexes_: List[int]
    node_ids_: set[int]
    load_indexes_: List[int]
    pmin_: float
    pmax_: float
    cross_nodes_indexes_: List[int]
    diff_: float

    def table_name(self):
        return "area2"

    def key_field(self):
        return "npa"

    def unit_name(self):
        return "территория"

    def __init__(
        self, index: int, id: int, name: str, pbase: float, territory_id: int = None
    ):
        super().__init__(index, id, name)
        self.pbase_ = pbase
        node_table = rastr.table("node")
        node_pn = node_table.column("pn")
        node_territory = node_table.column("npa")
        node_pnr = node_table.column("pnr")
        node_table.set_selection(f"(!sta)&({self.key_field()}={id})")
        node_data = node_table.data("ny")
        self.node_ids_ = set([node_id for node_id in node_data["ny"]])
        self.uncontrolled_load_ = 0.0
        self.lrc_load_ = 0.0
        self.node_indexes_ = []
        self.cross_nodes_indexes_ = []
        for node_id in self.node_ids_:
            node_table.set_selection(f"ny={node_id}")
            node_index = node_table.find_next(-1)
            self.node_indexes_.append(node_index)
            self.uncontrolled_load_ += node_pnr.z(node_index)
            self.lrc_load_ += node_pnr.z(node_index) - node_pn.z(node_index)
            if territory_id is not None:
                if territory_id == node_territory.z(node_index):
                    self.cross_nodes_indexes_.append(node_index)
        self.load_indexes_ = []
        self.pmin_ = 0.0
        self.pmax_ = 0.0
        self.pbase_by_loads_ = 0.0
        self.diff_ = 0

    def attach_load(self, load: Load, loads: LoadMap):
        if load.node_id_ in self.node_ids_:
            self.pbase_by_loads_ += load.pbase_
            if not load.const_:
                self.pmin_ += load.pmin_
                self.pmax_ += load.pmax_
                self.uncontrolled_load_ -= load.pbase_
                if load.id_ not in loads.load_keys_:
                    loads.loads_.append(load)
                    loads.load_keys_[load.id_] = len(loads.loads_) - 1
                self.load_indexes_.append(loads.load_keys_[load.id_])
            else:
                self.pmin_ += load.pbase_
                self.pmax_ += load.pbase_

    @staticmethod
    def dump_header():
        print(
            "Объединение                                   Нагр Ограничения         P0       Pсумнагр  Pнеконт   Pсхн   Pобщ"
        )

    def dump(self):
        print(
            f"{self.unit_name():12} {self.verbal_name()} "
            f"{len(self.load_indexes_):5} "
            f"[{self.pmin_:8.3f};{self.pmax_:8.3f}] "
            f"{self.pbase_:8.3f} "
            f"{self.pbase_by_loads_:8.3f} "
            f"{self.uncontrolled_load_:8.3f} "
            f"{self.lrc_load_:8.3f} "
            # f"{self.cross_load_:8.3f}"
        )

    def create_equation(self, loads: LoadMap):
        p = 0.0
        for load_index in self.load_indexes_:
            load = loads.loads_[load_index]
            p += load.pbase_
        print(f"{self.verbal_name()} {p} = {self.pbase_ - self.uncontrolled_load_}")
        return p

    def p_result(self):
        node_table = rastr.table("node")
        node_pn = node_table.column("pn")
        node_pnr = node_table.column("pnr")
        self.lrc_load_ = 0.0

        for node_index in self.node_indexes_:
            self.lrc_load_ += node_pnr.z(node_index) - node_pn.z(node_index)

        p_cross = 0.0
        for cross_node_index in self.cross_nodes_indexes_:
            p_cross += node_pnr.z(cross_node_index)

        p = rastr.table(self.table_name()).column("pn").z(self.index_)
        cross_coe = 1.0  # 1.2 * (1.0 - p_cross / p)
        self.diff_ = (self.pbase_ - p) * cross_coe
        return p


class Area(Territory):
    def table_name(self):
        return "area"

    def key_field(self):
        return "na"

    def unit_name(self):
        return "район"

    def __init__(
        self, index: int, id: int, name: str, pbase: float, territory_id: int = None
    ):
        super().__init__(index, id, name, pbase, territory_id)


model_file = "e:/downloads/4_1_11 Примеры/Пример_3_.os"
template_file = "d:/Documents/RastrWin3/SHABLON/poisk.os"
rastr = astra.Rastr()
rastr.load(astra.LoadCode.REPL, model_file, template_file)

territory_id = 52902

territory_table = rastr.table("area2")
territory_table.set_selection(f"npa={territory_id}")
territory_index = territory_table.find_next(-1)
if territory_index <= 0:
    raise Exception(f"Территория с номером {territory_id} не найдена")

areas = []

territory = Territory(
    territory_index,
    territory_id,
    territory_table.column("name").z(territory_index),
    territory_table.column("pn").z(territory_index),
)

areas.append(territory)

node_table = rastr.table("node")
node_area = node_table.column("na")
node_table.set_selection(f"npa={territory_id}")
index = node_table.find_next(-1)
area_ids = set()
while index >= 0:
    area_ids.add(node_area.z(index))
    index = node_table.find_next(index)

area_table = rastr.table("area")
area_id = area_table.column("na")
area_name = area_table.column("name")
area_pn = area_table.column("pn")


for id in area_ids:
    area_table.set_selection(f"na={id}")
    index = area_table.find_next(-1)
    if index >= 0:
        areas.append(
            Area(
                index,
                area_id.z(index),
                area_name.z(index),
                area_pn.z(index),
                territory_id,
            )
        )
    else:
        raise Exception(f"Район с номером {id} не найден")

load_table = rastr.table("Load")
load_table.set_selection("!sta")
load_data = zip(
    *rastr.table("Load").data("Num,Name,Node,Pmin,Pmax,P,Const,TPO").values()
)

load_map = LoadMap()

for (
    load_id,
    load_name,
    load_node,
    load_pmin,
    load_pmax,
    load_p,
    load_const,
    load_tpo,
    load_index,
) in load_data:
    load = Load(
        load_index,
        load_id,
        load_name,
        load_node,
        load_p,
        load_pmin,
        load_pmax,
        load_const,
    )
    for area in areas:
        area.attach_load(load, load_map)

Area.dump_header()
for area in areas:
    area.dump()

load_map.check_assert()

n_loads = len(load_map.loads_)
n_areas = len(areas)

# 1. Целевая функция: минимизируем sum((x_i - pbase_i)^2)
# Это эквивалентно 1/2 * x^T * P * x + q^T * x
# Где P = 2*I, q = -2*pbase
P = sparse.eye(n_loads, format="csc") * 2.0
q = np.array([-2.0 * load.pbase_ for load in load_map.loads_])

# 2. Матрица ограничений A (размер: n_areas + n_loads  x  n_loads)
# Первые n_areas строк — баланс по районам (sum(x_i) = Area_P_target)
# Остальные n_loads строк — границы каждой нагрузки (l_i <= x_i <= u_i)
A_rows = []
for area in areas:
    row = np.zeros(n_loads)
    for load_idx in area.load_indexes_:
        row[load_idx] = 1.0
    A_rows.append(row)

A = sparse.vstack(
    [sparse.csc_matrix(A_rows), sparse.eye(n_loads, format="csc")], format="csc"
)


tolerance = 0.2
territory.pbase_ = 1200

for iteration in range(1, 8):
    print(f"Итерация {iteration:4}")
    # 3. Векторы границ l и u
    # Для районов: равенство (l = u = P_target)
    # P_target = P_base_территории - неуправляемая_нагрузка
    area_targets = [a.pbase_ - a.uncontrolled_load_ + a.diff_ for a in areas]

    # Для нагрузок: [pmin, pmax]
    load_mins = [load.pmin_ for load in load_map.loads_]
    load_maxs = [load.pmax_ for load in load_map.loads_]

    lower = np.concatenate([area_targets, load_mins])
    upper = np.concatenate([area_targets, load_maxs])

    # 4. Решение
    prob = osqp.OSQP()
    prob.setup(P=P, q=q, A=A, l=lower, u=upper, verbose=False)
    res = prob.solve()

    # 5. Обработка результатов

    if res.info.status == "solved":
        load_p = load_table.column("P")
        for i, load in enumerate(load_map.loads_):
            # print(f"Нагрузка {load.id_}: было {load.pbase_:.3f} -> стало {res.x[i]:.3f}")
            load_p.set_z(load.index_, res.x[i])
        if rastr.rgm() == astra.ASTCode.OK:
            rastr.save("e:/temp/result.os", template_file)
            max_diff = 0.0
            for area, area_target in zip(areas, area_targets):
                p = area.p_result()
                diff = area.pbase_ - area.p_result()
                max_diff = max(abs(diff), max_diff)
                print(f"{area.verbal_name()} {p:8.3f} {diff:8.3f} {area.diff_:8.3f}")
            if max_diff < tolerance:
                print(f"Небаланс не превышает {tolerance}")
                break
        else:
            print("УР несбалансирован !")
            break
    else:
        print("Решение не найдено. Проверьте баланс мощностей и границы Pmin/Pmax.")
        break
